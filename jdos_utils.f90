!=============================================================================== 
! MODULE od_jdos  OptaDOS - Joint Density of States
! This is the module that contains all if the JDOS routines. It is used through 
! the global jdos_calculate subroutine
!-------------------------------------------------------------------------------
! NB. It should be possible to pass optioinal arguments to sub programs as
! optional argumnets without checking whether they are there or not. g95 will 
! allow this behaviour. gfotran will not.
!=============================================================================== 
module od_jdos_utils
  use od_algorithms, only : heap_sort,gaussian
  use od_constants,  only : bohr2ang, dp, H2eV
  use od_comms,      only : on_root, num_nodes, my_node_id, root_id,comms_slice,comms_bcast,& 
       &comms_send, comms_recv,comms_reduce
  use od_electronic, only : band_energy,band_gradient,efermi,efermi_castep,num_electrons, &
       &nbands, nspins, spin_polarised, elec_report_parameters, elec_read_band_gradient, &
       &elec_read_band_energy, electrons_per_state
  use od_io,         only : stdout,io_date,io_error,maxlen,seedname,filename_len, &
       &io_file_unit,io_time
  use od_cell,       only : real_lattice, recip_lattice, cell_volume, cell_find_MP_grid, & 
       &kpoint_r, kpoint_weight, nkpoints, kpoint_grid_dim, cell_calc_lattice, &
       &cell_report_parameters
  use od_dos_utils,        only : dos_utils_merge, doslin_sub_cell_corners, doslin
  use od_parameters

  implicit none


  !-------------------------------------------------------------------------------
  ! P U B L I C   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: jdos_adaptive(:,:)
  real(kind=dp), allocatable, public, save :: jdos_fixed(:,:)
  real(kind=dp), allocatable, public, save :: jdos_linear(:,:)

  integer,save :: jdos_nbins

  real(kind=dp), allocatable, public, save :: E(:)

  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S 
  public :: jdos_utils_calculate
  public :: write_jdos
  !-------------------------------------------------------------------------------

  real(kind=dp), save                   :: delta_bins ! Width of bins 
  logical :: calc_weighted_jdos
  integer, allocatable, save :: vb_max(:)
  !-------------------------------------------------------------------------------

contains

  !=============================================================================== 
  subroutine jdos_utils_calculate(matrix_weights, weighted_jdos) ! I've changed this
    !=============================================================================== 
    ! Main routine in dos module, drives the calculation of Density of states for
    ! both task : dos and also if it is required elsewhere.
    !=============================================================================== 
    implicit none
    integer :: ierr, idos, i, ik, is, ib
    real(kind=dp) :: time0, time1

    real(kind=dp),intent(out), allocatable, optional    :: weighted_jdos(:,:,:)  !I've added this
    real(kind=dp),intent(in), optional  :: matrix_weights(:,:,:,:)               !I've added this

    calc_weighted_jdos=.false.
    if(present(matrix_weights)) calc_weighted_jdos=.true.

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)')'+============================================================================+'
       write(stdout,'(1x,a78)')'+================== Joint Density Of States Calculation =====================+'
       write(stdout,'(1x,a78)')'+============================================================================+'
       write(stdout,*)
    endif

    if(calc_weighted_jdos.eqv..false.) then ! We are called just to provide dos.
       if(allocated(E)) then
          if(on_root) write(stdout,*) " Already calculated jdos, so returning..."
          return  ! The jdos has already been calculated previously so just return.       
       endif
    endif

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S 
    ! If we're using one of the more accurate roadening schemes we also need to read in the 
    ! band gradients too
    if(quad.or.linear.or.adaptive) call elec_read_band_gradient
    !-------------------------------------------------------------------------------


    allocate(vb_max(nspins), stat=ierr)
    if (ierr/=0) call io_error ("cannot allocate vb_max")


    ! For an insulator
    vb_max(:)=num_electrons(:)/electrons_per_state


    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   D O S 
    ! Now everything is set up, we can perform the dos accumulation in parellel
    time0=io_time()

    call setup_energy_scale

    if(on_root.and.(iprint>1)) write(stdout,*)

    if(fixed)then
       if(calc_weighted_jdos)then 
          call calculate_jdos(jdos_fixed, matrix_weights, weighted_jdos)
          call jdos_utils_merge(jdos_fixed,weighted_jdos)
       else
          call calculate_jdos(jdos_fixed)
          call jdos_utils_merge(jdos_fixed)
       endif


    endif
    if(adaptive)then
       if(calc_weighted_jdos)then 
          call calculate_jdos(jdos_adaptive, matrix_weights, weighted_jdos)
          call dos_utils_merge(jdos_adaptive, weighted_jdos)
       else
          call calculate_jdos(jdos_adaptive)
          call dos_utils_merge(jdos_adaptive)
       endif
    endif
    if(linear)then
       if(calc_weighted_jdos)then 
          call calculate_jdos(jdos_linear, matrix_weights, weighted_jdos)
          call dos_utils_merge(jdos_linear, weighted_jdos)
       else
          call calculate_jdos(jdos_linear)
          call dos_utils_merge(jdos_linear)
       endif
    endif


    if(quad) then
       call io_error("quadratic broadening not implemented") 
       !if(quad)    call merge_dos(dos_quad)
       !if(quad)    call merge_dos(intdos_quad)
    endif

    if(.not.on_root) then
       if(allocated(E)) deallocate(E, stat=ierr)
       if (ierr/=0) call io_error ("cannot deallocate  E")
    endif

    time1=io_time()
    if(on_root)  write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate dos  ',time1-time0,' (sec)'
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! F E R M I   E N E R G Y   A N A L Y S I S
    !time0=io_time()
    !write(stdout,*)
    !write(stdout,'(1x,a78)')  '+------------------------ Fermi Energy Analysis -----------------------------+'
    !write(stdout,'(1x,a1,a45,f8.4,a3,20x,a1)') "|","Fermi energy from CASTEP :",efermi_castep," eV","|"
    !write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
    !
    !if(compute_efermi) then
    ! if(fixed) then 
    !   write(stdout,'(1x,a78)') "| From Fixed broadening                                                      | "
    !   efermi_fixed= calc_efermi_from_intdos(intdos_fixed)
    !   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Fixed braodening) : ", Efermi_fixed,"eV","|"
    !   write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'

    ! endif
    ! if(adaptive) then
    !   write(stdout,'(1x,a78)') "| From Adaptive broadening                                                   | "  
    !   efermi_adaptive=calc_efermi_from_intdos(intdos_adaptive)
    !   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Adaptive braodening) : ", Efermi_adaptive,"eV","|"
    ! write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'

    ! endif
    ! if(linear) then
    !   write(stdout,'(1x,a78)') "| From Linear broadening                                                     | " 
    !   efermi_linear=calc_efermi_from_intdos(intdos_linear) 
    !   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Linear braodening) : ", Efermi_linear," eV","|"
    ! write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
    !
    ! endif
    !else
    ! write(stdout,'(1x,a78)')   "| No Fermi energies calculated                                               | "
    ! ! Use the derived Fermi energy shifts if calculated. It not use the fermi_energy supplied
    ! ! if not, use the CASTEP Fermi energy
    ! if(fermi_energy.ne.-990.0_dp)then
    !  efermi_fixed=fermi_energy
    !  efermi_adaptive=fermi_energy
    !  efermi_linear=fermi_energy
    ! else
    !  efermi_fixed=efermi_castep
    !  efermi_adaptive=efermi_castep
    !  efermi_linear=efermi_castep
    ! endif
    !endif

    ! NB If you have asked for more than one type of broadening
    ! If one of your options is linear then all will have the linear efermi
    ! If you have asked for adaptive, but nor linear, you will have the adaptive efermi
    !if(fixed) then
    ! efermi=efermi_fixed  
    !endif

    !if(adaptive) then
    ! efermi=efermi_adaptive
    !endif
    !
    !if(linear)then
    ! efermi=efermi_linear
    !endif
    !
    ! write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy used : ", efermi,"eV","|"
    ! E(:)=E(:)-efermi
    ! band_energy(:,:,:) = band_energy(:,:,:) - efermi
    !
    !write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
    !time1=io_time()
    !write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Fermi energies ',time1-time0,' (sec)'
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! B A N D   E N E R G Y   A N A L Y S I S
    ! Now for a bit of crosschecking  band energies
    ! These should all converge to the same number as the number of bins is increased
    !if(compute_band_energy) call compute_band_energies
    !-------------------------------------------------------------------------------


    if(dos_per_volume) then
       if(fixed) then
          jdos_fixed=jdos_fixed/cell_volume   
       endif
       if(adaptive) then
          jdos_adaptive=jdos_adaptive/cell_volume 
       endif
       if(linear) then
          jdos_linear=jdos_linear/cell_volume     
       endif

       ! if(quad) then
       !    dos_quad=dos_quad/cell_volume         
       !    intdos_quad=intdos_quad/cell_volume
       ! endif   
    endif

    if (on_root) then
       write(stdout,'(1x,a78)')    '+============================================================================+'
       write(stdout,'(1x,a78)')    '+============== Joint Density Of States Calculation End =====================+'
       write(stdout,'(1x,a78)')    '+============================================================================+'
       write(stdout,*)
    end if

  end subroutine jdos_utils_calculate


  !=============================================================================== 
  subroutine write_jdos(E,dos,dos_name)
    !=============================================================================== 
    ! This routine receives an energy scale, a density of states and a file name
    ! and writes out the DOS to disk
    !=============================================================================== 
    implicit none
    real(dp), intent(in) :: E(jdos_nbins)
    real(dp), intent(in) :: dos(jdos_nbins,nspins)
    character(len=*), intent(in) :: dos_name
    integer :: i, dos_file, ierr
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: dos_units, intdos_units


    dos_file=io_file_unit()
    open(unit=dos_file,file=trim(seedname)//'.j'//trim(dos_name)//'.dat',iostat=ierr)
    if(ierr.ne.0) call io_error(" ERROR: Cannot open output file in dos: write_dos")

    dos_units="(electrons per eV)" ; intdos_units="(electrons)"
    if(dos_per_volume) then
       dos_units="(electrons per eV/A^3)" 
       intdos_units="(electrons per A^3)"
    endif

    write(dos_file, *) "##############################################################################"
    write(dos_file,*) "#"
    write(dos_file, *) "#                  O p t a D O S   o u t p u t   f i l e "  
    write(dos_file, '(1x,a1)') "#"
    write(dos_file,*) "#    Denisty of States using ", trim(dos_name), " broadening"
    call io_date(cdate,ctime)
    write(dos_file,*)  '#  Generated on ',cdate,' at ',ctime
    write(dos_file,*) "# Column        Data"
    write(dos_file,*) "#    1        Energy (eV)"
    if(nspins>1) then
       write(dos_file,*) "#    1        Up-spin DOS ", trim(dos_units)
       write(dos_file,*) "#    2        D50own-spin DOS ", trim(dos_units)
       write(dos_file,*) "#    3        Up-spin Integrated DOS ", trim(intdos_units)
       write(dos_file,*) "#    4        Down-spin Integrated DOS ", trim(intdos_units)
    else
       write(dos_file,*) "#    2        DOS ", trim(dos_units)
       write(dos_file,*) "#    3        Integrated DOS ", trim(intdos_units)
    endif
    write(dos_file, '(1x,a1)') "#"
    write(dos_file, '(1x,a78)') "##############################################################################"


    if(nspins>1) then
       do i=1,jdos_nbins
          write(dos_file, *) E(i), dos(i,1), -dos(i,2)
       enddo
    else
       do i=1,jdos_nbins
          write(dos_file, *) E(i), dos(i,1)
       enddo
    endif
    close(dos_file)
  end subroutine write_jdos


  !=============================================================================== 
  subroutine setup_energy_scale
    !=============================================================================== 
    ! Sets up all broadening independent DOS concerns
    ! Calls the relevant dos calculator.
    !=============================================================================== 
    use od_dos_utils, only : dos_utils_calculate
    use od_parameters, only : compute_efermi
    use od_electronic, only : efermi_castep

    implicit none

    integer       :: idos,i,ierr
    real(kind=dp) :: max_band_energy

    if(compute_efermi)then
       call dos_utils_calculate
    else
       efermi=efermi_castep
    endif

    if(jdos_max_energy<0.0_dp) then ! we have to work it out ourselves
       max_band_energy=maxval(band_energy)
       call comms_reduce(max_band_energy,1,'MAX')
       call comms_bcast(max_band_energy,1)
       jdos_max_energy=efermi-max_band_energy

       if(on_root.and.(iprint>2)) then
          write(stdout,*)
          write(stdout,'(1x,a40,f11.3,a14)') 'max_band_energy (before correction) : ', max_band_energy, " <-- JDOS Grid"
       endif
    endif
    
    jdos_nbins=abs(ceiling(jdos_max_energy/jdos_spacing))
    jdos_max_energy=jdos_nbins*jdos_spacing
    
    allocate(E(1:jdos_nbins),stat=ierr)
    if (ierr/=0) call io_error ("cannot allocate E")

    delta_bins=jdos_max_energy/real(jdos_nbins-1,dp)
    do idos=1,jdos_nbins
       E(idos) = real(idos-1,dp)*delta_bins
    end do

  if(on_root.and.(iprint>2))then
     write(stdout,'(1x,a40,f11.5,a14)')  'efermi : ', efermi,  " <-- JDOS Grid"
       write(stdout,'(1x,a40,f11.5,a14)') 'jdos_max_energy : ', jdos_max_energy, " <-- JDOS Grid"
       write(stdout,'(1x,a40,i11,a14)') ' jdos_nbins : ', jdos_nbins, " <-- JDOS Grid"
       write(stdout,'(1x,a40,f11.3,a14)')' jdos_spacing : ', jdos_spacing, " <-- JDOS Grid"
       write(stdout,'(1x,a40,f11.3,a14)')' delta_bins : ', delta_bins, " <-- JDOS Grid"
    endif

  end subroutine setup_energy_scale



  !===============================================================================
  subroutine allocate_jdos(dos)
    !===============================================================================
    !===============================================================================
    implicit none

    real(kind=dp), allocatable  :: dos(:,:)
    real(kind=dp), allocatable  :: intdos(:,:)

    integer :: ierr 

    allocate(dos(jdos_nbins,nspins), stat=ierr)
    if(ierr/=0) call io_error("error in allocating dos")
    dos=0.0_dp

  end subroutine allocate_jdos


  !===============================================================================
  subroutine jdos_deallocate
    !===============================================================================
    !===============================================================================
    implicit none
    if(allocated(jdos_adaptive)) deallocate(jdos_adaptive)
    if(allocated(jdos_fixed))    deallocate(jdos_fixed)
    if(allocated(jdos_linear))   deallocate(jdos_linear) 
    if(allocated(E))             deallocate(E)
  end subroutine jdos_deallocate


  !===============================================================================
  subroutine calculate_jdos(jdos, matrix_weights, weighted_jdos)  ! I've changed this
    !===============================================================================

    !===============================================================================
    use od_comms, only : my_node_id, on_root
    use od_cell, only : num_kpoints_on_node
    implicit none

    integer :: i,ik,is,ib,idos,ierr,iorb,jb
    integer :: m,n,o,nn
    real(kind=dp) :: dos_temp, cuml, intdos_accum, width
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4)

    real(kind=dp),intent(inout),allocatable, optional    :: weighted_jdos(:,:,:)
    real(kind=dp),intent(in), optional  :: matrix_weights(:,:,:,:)

    real(kind=dp),intent(out),allocatable :: jdos(:,:)


    if(linear.or.adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:),dp)/2.0_dp
    if(adaptive) adaptive_smearing=adaptive_smearing*sum(step(:))/3
    if(fixed) width=fixed_smearing

    call allocate_jdos(jdos)
    if(calc_weighted_jdos) then
       allocate(weighted_jdos(jdos_nbins, nspins, 1))
       weighted_jdos=0.0_dp
    endif

    do ik=1,num_kpoints_on_node(my_node_id)
       if(iprint>1 .and. on_root) then
          if (mod(real(ik,dp),10.0_dp) == 0.0_dp) write(stdout,'(a40,i4,a3,i4,1x,a14,5x,a8)') &
               &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id),'on this node.',"<-- JDOS"
       endif
       do is=1,nspins
          do ib=1,vb_max(is)
             do jb=vb_max(is)+1,nbands
                if(linear.or.adaptive) grad(:) = real(band_gradient(jb,jb,:,ik,is)-band_gradient(ib,ib,:,ik,is),dp)
                if(linear) call doslin_sub_cell_corners(grad,step,band_energy(jb,is,ik)-band_energy(ib,is,ik),EV)
                if(adaptive) width = sqrt(dot_product(grad,grad))*adaptive_smearing

                ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
                ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
                if(finite_bin_correction.and.(width<delta_bins)) width = delta_bins

                do idos=1,jdos_nbins
                   ! The linear method has a special way to calculate the integrated dos
                   ! we have to take account for this here.
                   if(linear)then
                      dos_temp=doslin(EV(0),EV(1),EV(2),EV(3),EV(4),E(idos),cuml)
                   else
                      dos_temp=gaussian(band_energy(jb,is,ik)-band_energy(ib,is,ik),width,E(idos))&
                           &*electrons_per_state*kpoint_weight(ik)
                   endif

                   jdos(idos,is)=jdos(idos,is) + dos_temp  

                   ! this will become a loop over final index (polarisation)
                   ! Also need to remove kpoints weights.
                   if(calc_weighted_jdos) then
                      weighted_jdos(idos,is,1)=weighted_jdos(idos,is,1) + dos_temp*matrix_weights(ib,jb,ik,is)
                   end if
                      

                end do
             end do
          end do
       end do
    end do

  end subroutine calculate_jdos

  !=============================================================================== 
  subroutine jdos_utils_merge(jdos, weighted_jdos)
    !===============================================================================
    ! The DOS was calculated accross nodes. Now give them all back to root
    ! and free up the memeory on the slaves
    !------------------------------------------------------------------------------- 
    ! Arguments: dos          (in - slaves) (inout -  root)       : The DOS
    !            weighted_dos (in - slaves) (inout -  root) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !===============================================================================  
    use od_comms!,      only : on_root, comms_reduce
    use od_electronic, only : nspins
    use od_io,         only : io_error 

    implicit none
    integer :: idos,ierr
    real(kind=dp),intent(inout), allocatable, optional :: weighted_jdos(:,:,:) ! bins.spins, orbitals
    real(kind=dp),allocatable,intent(inout) :: jdos(:,:)

    call comms_reduce(jdos(1,1),nspins*jdos_nbins,"SUM")

    if(present(weighted_jdos))  call comms_reduce(weighted_jdos(1,1,1),nspins*jdos_nbins*1,"SUM")

    if(.not.on_root) then 
       if(allocated(jdos)) deallocate(jdos,stat=ierr)
       if (ierr/=0) call io_error (" ERROR : jdos : merge_jdos : cannot deallocate dos")
       if(present(weighted_jdos))  then
          if(allocated(weighted_jdos)) deallocate(weighted_jdos,stat=ierr)
          if (ierr/=0) call io_error (" ERROR : jdos : merge_jdos : cannot deallocate weighted_dos")
       end if
    endif
  end subroutine jdos_utils_merge



end module od_jdos_utils
