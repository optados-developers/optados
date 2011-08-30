!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 


module od_core

  use od_constants, only : dp
  implicit none
  private
  public :: core_calculate

  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_dos
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_dos_broadened



contains

  subroutine core_calculate
    use od_electronic, only  : elec_read_elnes_mat
    use od_dos_utils, only : dos_utils_calculate
    use od_comms, only : on_root
    use od_io, only : stdout
    use od_parameters, only : core_LAI_broadening, LAI_gaussian, LAI_lorentzian

    implicit none

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+============================ Core Loss Calculation =========================+'
       write(stdout,'(1x,a78)') '+============================================================================+' 
      write(stdout,*)
    endif


    ! read in the core matrix elements from disk
    call  elec_read_elnes_mat
    !    (elnes_mat(orb,nb,indx,nk,ns),indx=1,3)

    call core_prepare_matrix_elements

    call dos_utils_calculate(matrix_weights, weighted_dos)

    ! Lifetime and instrumental broadening
    if(core_LAI_broadening.eqv..true.) then 
       allocate(weighted_dos_broadened(size(weighted_dos,1),size(weighted_dos,2),size(weighted_dos,3)))
       weighted_dos_broadened=0.0_dp
       if(LAI_gaussian) call core_gaussian 
       if(LAI_lorentzian) call core_lorentzian 
    endif

    if (on_root) then
       call write_core
    endif

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+========================== Core Loss Calculation End =======================+'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif


  end subroutine core_calculate

  ! Private routines

  subroutine core_prepare_matrix_elements
    use od_electronic, only  : elnes_mat,  elnes_mwab, nbands, nspins,num_electrons, electrons_per_state
    use od_comms, only : my_node_id
    use od_cell, only : num_kpoints_on_node
    use od_parameters, only : core_geom, core_qdir ! needs changing to elnes
    use od_io, only : io_error

    real(kind=dp), dimension(3) :: qdir 
    real(kind=dp) :: q_weight
    integer :: N, N_spin, n_eigen, orb, ierr
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp) :: g

    num_occ = 0.0_dp
    do N_spin=1,nspins
       num_occ(N_spin) = num_electrons(N_spin)
    enddo

    if (electrons_per_state==2) then 
       num_occ(1) = num_occ(1)/2.0_dp
    endif

    allocate(matrix_weights(elnes_mwab%norbitals,elnes_mwab%nbands,num_kpoints_on_node(my_node_id),nspins),stat=ierr)
    if(ierr/=0) call io_error('Error: core_prepare_matrix_elements - allocation failed for matrix_weights')
    matrix_weights=0.0_dp

    if (index(core_geom,'polar')>0) then 
       qdir=core_qdir            
       q_weight=((qdir(1)**2.0_dp)+(qdir(2)**2.0_dp)+(qdir(3)**2.0_dp))**0.5_dp
       if(q_weight<0.001_dp)&
            call io_error("Error: core_prepare_matrix_elements.  please check core_qdir, norm close to zero")
    end if

    do N=1,num_kpoints_on_node(my_node_id)                      ! Loop over kpoints
       do N_spin=1,nspins                                    ! Loop over spins
          do n_eigen=(nint(num_occ(N_spin))+1),nbands    ! Loop over unoccupied states
             do orb=1,elnes_mwab%norbitals
                if(index(core_geom,'polar')>0) then 
                   g = (((qdir(1)*elnes_mat(orb,n_eigen,1,N,N_spin))+ &
                        (qdir(2)*elnes_mat(orb,n_eigen,2,N,N_spin))+ &
                        (qdir(3)*elnes_mat(orb,n_eigen,3,N,N_spin)))/q_weight)
                   matrix_weights(orb,n_eigen,N,N_spin) = real(g*conjg(g),dp)
                else
                   matrix_weights(orb,n_eigen,N,N_spin) = ( &
                        elnes_mat(orb,n_eigen,1,N,N_spin)*conjg(elnes_mat(orb,n_eigen,1,N,N_spin)) + &
                        elnes_mat(orb,n_eigen,2,N,N_spin)*conjg(elnes_mat(orb,n_eigen,2,N,N_spin)) + &
                        elnes_mat(orb,n_eigen,3,N,N_spin)*conjg(elnes_mat(orb,n_eigen,3,N,N_spin)) ) / 3.0_dp
                   !                matrix_weights(orb,n_eigen,N,N_spin) = real(g*conjg(g),dp)
                   !           matrix_weights(n_eigen,n_eigen2,N,N_spin) = 1.0_dp  ! 
                end if
             end do
          end do
       end do
    end do

  end subroutine core_prepare_matrix_elements

  subroutine write_core
    !***************************************************************
    ! This subroutine writes out the Core loss function

    use od_constants, only : bohr2ang, periodic_table_name
    use od_parameters, only : dos_nbins, core_LAI_broadening, LAI_gaussian, LAI_gaussian_width, &
         LAI_lorentzian, LAI_lorentzian_scale,  LAI_lorentzian_width,  LAI_lorentzian_offset, output_format
    use od_electronic, only: elnes_mwab,elnes_orbital
    use od_io, only : seedname, io_file_unit,io_error
    use od_dos_utils, only : E
    use od_cell, only : num_species, atoms_symbol
    use xmgrace_utils

    integer :: N
    real(kind=dp) ::dE,min_x,min_y,max_x,max_y, range
    integer :: core_unit,orb,ierr,loop,loop2,counter,num_edge,num_sites
    character(len=20) :: temp
    character(len=40) :: temp2
    character(len=10), allocatable :: elnes_symbol(:)
    integer, allocatable :: edge_shell(:), edge_am(:), edge_num_am(:), edge_list(:,:)
    integer, allocatable :: edge_species(:),edge_rank_in_species(:)
    integer, allocatable :: ion_species(:),ion_num_in_species(:)
    character(len=40), allocatable :: edge_name(:)
    real(kind=dp), allocatable :: dos_temp(:), dos_temp2(:)
    logical :: found

    allocate(dos_temp(dos_nbins),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of dos_temp failed')
    allocate(dos_temp2(dos_nbins),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of dos_temp2 failed')
    allocate(elnes_symbol(maxval(elnes_orbital%species_no(:))),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of elnes_symbol failed')

    dE=E(2)-E(1)

    counter=1
    do loop2=1,109
       do loop=1,num_species
          if(atoms_symbol(loop)==periodic_table_name(loop2)) then
             elnes_symbol(counter)=periodic_table_name(loop2)
             counter=counter+1
             !check atom count here
          end if
       end do
    end do

    ! Open the output file
    core_unit = io_file_unit()
    open(unit=core_unit,action='write',file=trim(seedname)//'_core.dat')

    ! Write into the output file
    write(core_unit,*)'#*********************************************'
    write(core_unit,*)'#            Core loss function               '
    write(core_unit,*)'#*********************************************'
    write(core_unit,*)'#'
    if(core_LAI_broadening) then 
       if(LAI_gaussian) write(core_unit,*)'# Gaussian broadening: FWHM', LAI_gaussian_width
       if(LAI_lorentzian) then 
          write(core_unit,*)'# Lorentzian broadening included'
          write(core_unit,*)'# Lorentzian scale ', LAI_lorentzian_scale
          write(core_unit,*)'# Lorentzian offset ', LAI_lorentzian_offset
          write(core_unit,*)'# Lorentzian width ', LAI_lorentzian_width
       end if
    end if
    write(core_unit,*)'#'

    weighted_dos=weighted_dos*bohr2ang**2  ! Converts units, note I don't have to worry about this in optics.f90 at electronic does it 
    if (core_LAI_broadening) weighted_dos_broadened=weighted_dos_broadened*bohr2ang**2

    do orb=1,elnes_mwab%norbitals   ! writing out doesn't include spin at the moment. 

       write(temp,'(a2,i1,a17)') 'n=',elnes_orbital%shell(orb),' ang= '//trim(elnes_orbital%am_channel_name(orb))

       write(temp2,'(a5,a2,1x,i0,a28)') 'Ion: ',trim(elnes_symbol(elnes_orbital%species_no(orb))),&
	& elnes_orbital%rank_in_species(orb),' State: '//trim(temp)
       write(core_unit,*) '# ',trim(temp2)


       do N=1,dos_nbins
          if (core_LAI_broadening) then 
             write(core_unit,*)E(N),weighted_dos(N,1,orb),weighted_dos_broadened(N,1,orb)
          else 
             write(core_unit,*)E(N),weighted_dos(N,1,orb)
          end if
       end do
       write(core_unit,*)''
    end do

    close(core_unit)

    allocate(ion_species(elnes_mwab%norbitals))
    allocate(ion_num_in_species(elnes_mwab%norbitals))
    ion_species=0;ion_num_in_species=0
    !Find out what ions we have
    do loop=1,elnes_mwab%norbitals
       if(loop==1) then
          ion_species(1)=elnes_orbital%species_no(1)
          ion_num_in_species(1)=elnes_orbital%rank_in_species(1)
          counter=1
       else
          found=.false.
          do loop2=1,counter
             if(elnes_orbital%species_no(loop)==ion_species(loop2).and.&
	& elnes_orbital%rank_in_species(loop)==ion_num_in_species(loop2)) then
                found=.true.
             end if
          enddo
          if(.not.found) then
             counter=counter+1
             ion_species(counter)=elnes_orbital%species_no(loop)
             ion_num_in_species(counter)=elnes_orbital%rank_in_species(loop)
          endif
       endif
    end do
    num_sites=counter

    ! We allocate these arrays as the max possible size, and just fill in the bits we need
    allocate(edge_species(elnes_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_species failed')
    allocate(edge_rank_in_species(elnes_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_rank_in_species failed')
    allocate(edge_shell(elnes_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_shell failed')
    allocate(edge_am(elnes_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_am failed')
    allocate(edge_num_am(elnes_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_num_am failed')
    allocate(edge_list(elnes_mwab%norbitals,7),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_list failed')
    edge_species=0;edge_rank_in_species=0;edge_shell=0;edge_am=0;edge_num_am=0;edge_list=0

    counter=1
    ! Find out how many edges
    do loop=1,elnes_mwab%norbitals
       if(loop==1) then
          edge_species(counter)=elnes_orbital%species_no(loop)
          edge_rank_in_species(counter)=elnes_orbital%rank_in_species(loop)
          edge_shell(counter)=elnes_orbital%shell(loop)
          edge_am(counter)=elnes_orbital%am_channel(loop)
          edge_num_am(counter)=edge_num_am(counter)+1
          edge_list(counter,edge_num_am(counter))=loop
       else
          ! else check if we have this am state
          found=.false.
          do loop2=1,counter
             if(edge_species(loop2)==elnes_orbital%species_no(loop).and.&
	& edge_rank_in_species(loop2)==elnes_orbital%rank_in_species(loop).and.&
                  edge_shell(loop2)==elnes_orbital%shell(loop).and.edge_am(loop2)==elnes_orbital%am_channel(loop)) then
                edge_num_am(counter)=edge_num_am(counter)+1
                edge_list(counter,edge_num_am(counter))=loop
                found=.true.
             end if
          end do
          if(.not.found) then
             counter=counter+1
             edge_species(counter)=elnes_orbital%species_no(loop)
             edge_rank_in_species(counter)=elnes_orbital%rank_in_species(loop)
             edge_shell(counter)=elnes_orbital%shell(loop)
             edge_am(counter)=elnes_orbital%am_channel(loop)
             edge_num_am(counter)=edge_num_am(counter)+1
             edge_list(counter,edge_num_am(counter))=loop
          end if
       end if
    end do
    num_edge=counter
    !
    allocate(edge_name(num_edge),stat=ierr)
    if(ierr/=0) call io_error('Error: core_write - allocation of edge_name failed')
    ! fill in edge name
    do loop=1,num_edge
       if(edge_shell(loop)==1) then
          temp='K1'
       elseif(edge_shell(loop)==2) then
          if(edge_am(loop)==0) then
             temp='L1'
          elseif(edge_am(loop)==1) then
             temp='L2,3'
          endif
       elseif(edge_shell(loop)==3) then
          if(edge_am(loop)==0) then
             temp='M1'
          elseif(edge_am(loop)==1) then
             temp='M2,3'
          elseif(edge_am(loop)==2) then
             temp='M4,5'
          endif
       elseif(edge_shell(loop)==4) then
          if(edge_am(loop)==0) then
             temp='N1'
          elseif(edge_am(loop)==1) then
             temp='N2,3'
          elseif(edge_am(loop)==2) then
             temp='N4,5'
          elseif(edge_am(loop)==3) then
             temp='N6,7'
          endif
       elseif(edge_shell(loop)==5) then
          if(edge_am(loop)==0) then
             temp='O1'
          elseif(edge_am(loop)==1) then
             temp='O2,3'
          elseif(edge_am(loop)==2) then  ! after this point I think we've drifted beyond what is physical!
             temp='O4,5'
          elseif(edge_am(loop)==3) then
             temp='O6,7'
          endif
       elseif(edge_shell(loop)==6) then
          if(edge_am(loop)==0) then
             temp='P1'
          elseif(edge_am(loop)==1) then
             temp='P2,3'
          elseif(edge_am(loop)==2) then
             temp='P4,5'
          elseif(edge_am(loop)==3) then
             temp='P6,7'
          endif
       endif
       
       write(edge_name(loop),'(a2,1x,i0,1x,a5)') trim(elnes_symbol(edge_species(loop))),edge_rank_in_species(loop),trim(temp)
    end do


    ! Now we know how many edges we have we can write them to a file


    ! Open the output file
    core_unit = io_file_unit()
    open(unit=core_unit,action='write',file=trim(seedname)//'_core_edge.dat')

    ! Write into the output file
    write(core_unit,*)'#*********************************************'
    write(core_unit,*)'#            Core loss function               '
    write(core_unit,*)'#*********************************************'
    write(core_unit,*)'#'
    if(core_LAI_broadening) then 
       if(LAI_gaussian) write(core_unit,*)'# Gaussian broadening: FWHM', LAI_gaussian_width
       if(LAI_lorentzian) then 
          write(core_unit,*)'# Lorentzian broadening included'
          write(core_unit,*)'# Lorentzian scale ', LAI_lorentzian_scale
          write(core_unit,*)'# Lorentzian offset ', LAI_lorentzian_offset
          write(core_unit,*)'# Lorentzian width ', LAI_lorentzian_width
       end if
    end if
    write(core_unit,*)'#'

    do loop=1,num_edge
       write(core_unit,*) '# ',trim(edge_name(loop))

       dos_temp=0.0_dp;dos_temp2=0.0_dp
       do loop2=1,edge_num_am(loop)
          dos_temp=dos_temp+weighted_dos(:,1,edge_list(loop,loop2))/real(edge_num_am(loop),dp)
          if (core_LAI_broadening) then 
             dos_temp2=dos_temp2+weighted_dos_broadened(:,1,edge_list(loop,loop2))/real(edge_num_am(loop),dp)
          end if
       end do

       do N=1,dos_nbins
          if (core_LAI_broadening) then 
             write(core_unit,*)E(N),dos_temp(N),dos_temp2(N)
          else 
             write(core_unit,*)E(N),dos_temp(N)
          end if
       end do
       write(core_unit,*)''
    end do

    close(core_unit)



    if(num_sites==1) then ! if only one site we write out plot script files

       if(trim(output_format)=="xmgrace") then
       
          core_unit=io_file_unit()
          open(unit=core_unit,file=trim(seedname)//'_'//'core_edge'//'.agr',iostat=ierr)
          if(ierr.ne.0) call io_error(" ERROR: Cannot open xmgrace batch file in core: write_core_xmgrace")

          min_x=minval(E)
          max_x=maxval(E)

          min_y=minval(weighted_dos)
          max_y=maxval(weighted_dos)


          ! For aesthetic reasons we make the axis range 1% larger than the data range
          range=abs(max_y-min_y)
          max_y=max_y+0.01_dp*range
          min_y=0.0_dp!min_y-0.01_dp*range

          call  xmgu_setup(core_unit)
          call  xmgu_legend(core_unit)
          call  xmgu_title(core_unit, min_x, max_x, min_y, max_y, 'Core-loss Spectrum')

          call  xmgu_axis(core_unit,"y",'Units')
          call  xmgu_axis(core_unit,"x",'Energy (eV)')

          do loop=1,num_edge
             dos_temp=0.0_dp;dos_temp2=0.0_dp
             do loop2=1,edge_num_am(loop)
                dos_temp=dos_temp+weighted_dos(:,1,edge_list(loop,loop2))/real(edge_num_am(loop),dp)
             end do

             call xmgu_data_header(core_unit,loop,loop,trim(edge_name(loop)))
             call xmgu_data(core_unit,loop,E(:),dos_temp)
          end do
       endif

    end if




  end subroutine write_core

  subroutine core_gaussian 
    !**************************************************************
    ! This subroutine adds in instrumental (Gaussian) broadening  

    use od_constants, only : pi, dp, bohr2ang 
    use od_parameters, only : LAI_gaussian_width, dos_nbins
    use od_dos_utils, only : E
    use od_electronic, only : nspins, elnes_mwab, band_energy

    integer :: N, N_spin, N_energy, N_energy2 
    real(kind=dp) :: G_width, g, dE

    G_width = LAI_gaussian_width          ! FWHM of Gaussian 
    dE = E(2)-E(1)

    do N=1,elnes_mwab%norbitals         ! Loop over orbitals 
       do N_spin=1,nspins               ! Loop over spins 
          do N_energy=1,dos_nbins       ! Loop over energy 
             do N_energy2=1,dos_nbins   ! Turn each energy value into a function 
                g = (((4.0_dp*log(2.0_dp))/pi)**(0.5_dp))*(1/G_width)*exp(-4.0_dp*(log(2.0_dp))*&
                     (((E(N_energy2)-E(N_energy))/G_width)**2.0_dp))  ! Gaussian 
                weighted_dos_broadened(N_energy2,N_spin,N)=weighted_dos_broadened(N_energy2,N_spin,N) &
                     + (g*weighted_dos(N_energy,N_spin,N)*dE)
             end do
          end do                        ! End loop over energy 
       end do                           ! End loop over spins  
    end do                              ! End loop over orbitals  

  end subroutine core_gaussian

  subroutine core_lorentzian 
    !**************************************************************
    ! This subroutine adds in life-time (Lorentzian) broadening  

    use od_constants, only : pi, dp
    use od_parameters, only : LAI_lorentzian_width, LAI_lorentzian_scale, LAI_lorentzian_offset, &
         LAI_gaussian_width, dos_nbins, LAI_gaussian, compute_efermi, adaptive, linear, fixed
    use od_dos_utils, only : E
    use od_electronic, only : nspins, elnes_mwab
    use od_dos_utils, only : efermi_fixed, efermi_adaptive, efermi_linear

    integer :: N, N_spin, N_energy, N_energy2 
    real(kind=dp) :: L_width, l, dE, e_fermi
    real(kind=dp),allocatable, dimension(:,:,:) :: weighted_dos_temp

    dE = E(2)-E(1)   
    allocate(weighted_dos_temp(size(weighted_dos,1),size(weighted_dos,2),size(weighted_dos,3)))
    weighted_dos_temp=0.0_dp

    if(LAI_gaussian) then 
       weighted_dos_temp = weighted_dos_broadened
       weighted_dos_broadened = 0.0_dp
    else 
       weighted_dos_temp = weighted_dos
    end if

    if(compute_efermi) then !this is assuming that I have already run the DOS!!!!  
       if(adaptive) e_fermi = efermi_adaptive
       if(linear) e_fermi = efermi_linear
       if(fixed)  e_fermi = efermi_fixed
       !   else call io_error ("OPTICS: No Ef set")
    endif

    do N=1,elnes_mwab%norbitals         ! Loop over orbitals 
       do N_spin=1,nspins               ! Loop over spins 
          do N_energy=1,dos_nbins       ! Loop over energy 
             if(E(N_energy).ge.(LAI_lorentzian_offset+e_fermi)) then 
                L_width = 0.5_dp*LAI_lorentzian_width & ! HWHW of Lorentzian 
                     + (E(N_energy)*LAI_lorentzian_scale)
             else 
                L_width = 0.5_dp*LAI_lorentzian_width 
             end if
             if (L_width.eq.0.0_dp) then ! if there is no broadening ie width=0 and below offset
                weighted_dos_broadened(N_energy,N_spin,N)=weighted_dos_temp(N_energy,N_spin,N) 
             else                        ! there is broadening 
                do N_energy2=1,dos_nbins ! Turn each energy value into a function 
                   l = weighted_dos_temp(N_energy,N_spin,N)*L_width/(pi*(((E(N_energy2)-E(N_energy))**2)+(L_width**2)))  ! Lorentzian 
                   weighted_dos_broadened(N_energy2,N_spin,N)=weighted_dos_broadened(N_energy2,N_spin,N)+(l*dE) 
                end do
             end if
          end do                        ! End look over energy 
       end do                           ! End loop over spins  
    end do                              ! End loop over orbitals  

  end subroutine core_lorentzian


!!$
!!$  !=============================================================================== 
!!$  subroutine write_core_gnuplot(label,E,column1,column2,column3)
!!$    !=============================================================================== 
!!$    use od_io,         only : io_file_unit,io_error,seedname 
!!$    implicit none 
!!$
!!$    type(graph_labels),intent(in) :: label
!!$
!!$    real(dp),  intent(in) :: E(:)
!!$    real(dp),  intent(in)  :: column1(:)
!!$    real(dp),  optional, intent(in) :: column2(:)
!!$    real(dp),  optional, intent(in) :: column3(:)
!!$
!!$    integer :: gnu_unit,ierr
!!$
!!$    gnu_unit=io_file_unit()
!!$    open(unit=gnu_unit,file=trim(seedname)//'_'//trim(label%name)//'.gnu',iostat=ierr)
!!$    if(ierr.ne.0) call io_error(" ERROR: Cannot open gnuplot batch file in optics: write_optics_gnupot")
!!$
!!$    gnu_unit = io_file_unit()
!!$    open(unit=gnu_unit,action='write',file=trim(seedname)//'_'//trim(label%name)//'.gnu')
!!$    write(gnu_unit,*) 'set xlabel ','"'//trim(label%x_label)//'"'
!!$    write(gnu_unit,*) 'set ylabel ','"'//trim(label%y_label)//'"'
!!$    write(gnu_unit,*) 'set title ','"'//trim(label%title)//'"'
!!$    if(present(column3)) then
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:3 t ','"'//trim(label%legend_b)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:4 t ','"'//trim(label%legend_c)//'"',' w l'
!!$    elseif(present(column2)) then
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:3 t ','"'//trim(label%legend_b)//'"',' w l'
!!$    else
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l'
!!$    endif
!!$    close(gnu_unit)
!!$
!!$  end subroutine write_optics_gnuplot
!!$

end module od_core
