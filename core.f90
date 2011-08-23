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

    use od_constants, only : bohr2ang
    use od_parameters, only : dos_nbins, core_LAI_broadening, LAI_gaussian, LAI_gaussian_width, &
         LAI_lorentzian, LAI_lorentzian_scale,  LAI_lorentzian_width,  LAI_lorentzian_offset
    use od_electronic, only: elnes_mwab
    use od_io, only : seedname, io_file_unit
    use od_dos_utils, only : E

    integer :: N
    real(kind=dp) ::dE
    integer :: core_unit,orb

    dE=E(2)-E(1)

    ! Open the output file
    core_unit = io_file_unit()
    open(unit=core_unit,action='write',file=trim(seedname)//'.core')

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
       DO N=1,dos_nbins
          if (core_LAI_broadening) then 
             write(core_unit,*)E(N),weighted_dos(N,1,orb),weighted_dos_broadened(N,1,orb)
          else 
             write(core_unit,*)E(N),weighted_dos(N,1,orb)
          end if
       end do
       write(core_unit,*)''
       write(core_unit,*)''
    end do

    close(core_unit)


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
          end do                        ! End look over energy 
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

end module od_core
