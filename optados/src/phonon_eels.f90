
!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!
! This file is part of OptaDOS
!
! OptaDOS - For obtaining electronic structure properties based on
!             integrations over the Brillouin zone
! Copyright (C) 2011  Andrew J. Morris,  R. J. Nicholls, C. J. Pickard
!                         and J. R. Yates
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!=========================================================================!
! MODULE od_phonon_eels
! This module contains routines for calculating a phonon EELS spectrum
!-------------------------------------------------------------------------------
module od_phonon_eels

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  use od_constants, only: dp
  use od_comms, only: on_root

  implicit none

  real(kind=dp), public, allocatable, save  :: temp(:, :, :)

  logical, public, save :: aloof_scattering
  logical, public, save :: impact_scattering
  logical, public, save :: semiclassical_aloof
  logical, public, save :: dipole_aloof

  !-------------------------------------------------------------------------!

  private

  public :: phonon_eels_calculate

  real(kind=dp), allocatable, save  :: adf(:, :) !(iatom,1:6)
  real(kind=dp), allocatable, save  :: debye_waller(:) !(iatom)
  complex(kind=dp), allocatable, save  :: phonon_eigenvectors(:, :, :, :) ! iqpoint, ieigenvalues, iatom, i=1,3)
  real(kind=dp), allocatable, save  :: phonon_eigenvalues(:, :)
  real(kind=dp), allocatable, save :: qpoint_positions(:, :)
    real(dp), allocatable :: atomic_positions(:, :) ! one day this ought to go in the right module, this includes mass as the 4th coordinate 
  
  real(kind=dp), allocatable, save  :: lf_dielectric_tensor(:, :, :) ! low frequency dielectric tensor

  real(kind=dp), save  :: phonon_lattice(1:3, 1:3) ! Maybe this is worth keeping

  real(kind=dp), save  :: inf_dielectric_tensor(1:3, 1:3) !  dielectric tensor  ^ infinity
  real(kind=dp), allocatable, save :: born_effective(:, :, :)
  real(kind=dp), allocatable, save :: mode_osc(:, :, :)
  
  integer :: num_ions ! found in phonon file
  integer :: num_eigenvalues
  integer :: num_qpoints

  integer, save :: nbins ! need to define this properly, should this go here?
  real(kind=dp), allocatable, save :: freq_scale(:) ! need to check what this is, should this go here?
  real(kind=dp), save :: freq ! need to move and think about energy as well

  
contains

  subroutine phonon_eels_calculate
    use od_io, only: stdout, io_error
    use od_electronic, only: elec_read_optical_mat, elec_read_band_energy
    use od_parameters, only: iprint, phonon_eels_task, phonon_eels_aloof_method

    implicit none

    integer :: ierr
    integer :: ibin ! need to think about where this should go
    real(kind=dp) :: broadening  ! need to think about where this should go
    real(kind=dp) :: kx  ! need to think about where this should go, what is it!
    real(kind=dp) :: ky  ! need to think about where this should go, what is it!
    real(kind=dp) :: impact_parameter ! This probably needs to be moved
    integer :: N
 
    
    aloof_scattering = .false.
    impact_scattering = .false.
    dipole_aloof = .false.
    semiclassical_aloof = .false.

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                   P h o n o n   E E L S  Calculation                       +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)')
    end if

    ! read what task we're supposed to be doing.
    call phonon_eels_read_task ! Completed

    ! Need to think about the frequency scale freq(1,nomega)
    ! O to maximum mode frequency+10%

  !  if (impact_scattering) then
      !  call phonon_eels_get_thermal_noise() ! PE_read_adf && PE_make_dewall or PE_read_debwall
                                              !  Completed       DUMMY             DUMMY
      !    call elec_read_band_energy()   ! Completed - shouldn't need this 
   !   call phonon_eels_read_phonon_file() ! Completed
   !   call phonon_eels_read_chge_trans()  ! DUMMY Subroutine
   !  endif

    if (aloof_scattering) then
      call phonon_eels_read_aloof_method   ! Completed
      if (dipole_aloof) then
        call io_error(' phonon_eels_calculate : Dipole Aloof Not yet implememnted.')
      elseif (semiclassical_aloof) then
        call phonon_eels_read_phonon_file() ! Completed

        freq=1000.0_dp ! This should be read in 
        nbins=50       ! This should be read in 

        allocate(freq_scale(1:nbins), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate lf_dielectric_tensor")
        freq_scale=0.0_dp

        call phonon_eels_set_frequency_scale()
        call phonon_eels_read_inf_dielectric_tensor()   

        allocate(born_effective(3,3,1:num_ions), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate lf_dielectric_tensor")
        born_effective=0.0_dp
        
        call phonon_eels_read_born_eff_charges()

        allocate(mode_osc(3,3,1:num_eigenvalues), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate lf_dielectric_tensor")
        mode_osc=0.0_dp

        call phonon_eels_calculate_mode_osc()
        broadening=0.01_dp ! probably needs to be in the input file

        ! Populate the lf_dielectric_tensor
        allocate(lf_dielectric_tensor(1:3, 1:3, 1:nbins), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate lf_dielectric_tensor")
        lf_dielectric_tensor=0.0_dp

        call phonon_eels_calculate_lf_dielectric_tensor(lf_dielectric_tensor, inf_dielectric_tensor, broadening)
         impact_parameter = 30.0_dp
 !       kx=0.01_dp
 !       ky=0.01_dp

   if (on_root) then
      write (stdout, '(1x,a1,a26,1x,f10.2,10x,a1)') "|", "Impact parameter (nm) : ", impact_parameter, "|"
      write (stdout, '(1x,a1,a26,1x,f10.2,10x,a1)') "|", 'Broadening_parameter :', broadening, "|"
      write (stdout, '(1x,a78)') 
      write (stdout, '(1x,a78)')
    end if       
        
  !       write(*,*) polarizability(ibin, kx, ky)

     else
        call io_error(' phonon_eels_calculate : Unknown aloof method.')
      endif
    endif


  end subroutine phonon_eels_calculate

  subroutine phonon_eels_set_frequency_scale ! Need to tidy this up to include freq/nbands.neq.integer etc
    use od_io, only: stdout
    implicit none
    real(kind=dp) :: dE
    integer :: N 
    
    dE = freq/nbins
     do N=1,nbins
        freq_scale(N)=N*dE
     end do 

  end subroutine  phonon_eels_set_frequency_scale

  subroutine phonon_eels_read_inf_dielectric_tensor
    use od_io, only: stdout, io_file_unit, io_error
    use od_io, only: seedname
    use od_parameters, only: iprint
    implicit none  
    integer :: phonon_in_unit, ierr, i
   
    phonon_in_unit = io_file_unit()
    open (unit=phonon_in_unit, file=trim(seedname)//".inf_dielectric_tensor", form='formatted', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open .inf_dielectric_tensor file in phonon_eels_read_inf_dielectric_tendor") ! this doesn't work....

    read (phonon_in_unit, *) (inf_dielectric_tensor(1, i), i=1, 3)
    read (phonon_in_unit, *) (inf_dielectric_tensor(2, i), i=1, 3)
    read (phonon_in_unit, *) (inf_dielectric_tensor(3, i), i=1, 3)
 !    write (stdout, *) "Read phonon_eels_read_dielectric_perm"   
    
  endsubroutine phonon_eels_read_inf_dielectric_tensor 

  subroutine phonon_eels_read_born_eff_charges
    use od_io, only: stdout, io_file_unit, io_error, maxlen
    use od_io, only: seedname
    use od_parameters, only: iprint
    implicit none  
    integer :: phonon_in_unit, ierr, i, N
    character(maxlen) :: dummy 
   
    phonon_in_unit = io_file_unit()
    open (unit=phonon_in_unit, file=trim(seedname)//".born_eff_charges", form='formatted', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open .born_eff_charges file in phonon_eels_read_born_eff_charges") ! this doesn't work....

    do N=1,num_ions
       read (phonon_in_unit, *) dummy, dummy, (born_effective(1, i, N), i=1, 3)
       read (phonon_in_unit, *) (born_effective(2, i, N), i=1, 3)
       read (phonon_in_unit, *) (born_effective(3, i, N), i=1, 3)
    end do

 !    write (stdout, *) "Read phonon_eels_read_dielectric_perm"   
    
  endsubroutine phonon_eels_read_born_eff_charges

  subroutine phonon_eels_calculate_mode_osc()
    implicit none
    real(kind=dp), allocatable :: mode_osc_1(:, :)
    real(kind=dp), allocatable :: mode_eff_charge(:, :)
    integer :: N, M, a, j
    real(kind=dp), allocatable :: norm(:)
    real(kind=dp) :: one_debye_to_eV ! might need to move this further up, is it used anywhere else?

    one_debye_to_eV = 0.2081943_dp
    
    allocate(mode_osc_1(3,1:num_eigenvalues))
    mode_osc_1=0.0_dp
    allocate(mode_eff_charge(3,1:num_eigenvalues))
    mode_eff_charge=0.0_dp 
    allocate(norm(1:num_eigenvalues))
    norm = 0.0_dp  
    
!    M=3
     do M=1,num_eigenvalues
       do N=1,num_ions
          do a=1,3
             do j=1,3
              !  write(*,*) a, N, j, mode_osc_1(a,1)
             !   mode_osc_1(a,M) = mode_osc_1(a,M) + born_effective(a,j,N)*real(phonon_eigenvectors(1,M,N,j))
                mode_osc_1(a,M) = mode_osc_1(a,M) + &
                born_effective(a,j,N)*(real(phonon_eigenvectors(1,M,N,j),dp))/sqrt(atomic_positions(N,4))
                 !                  write(*,*) a, N, j, mode_osc_1(a,1)
             end do
 !!         norm(M)=norm(M)+(real(phonon_eigenvectors(1,M,N,a),dp)*(real(phonon_eigenvectors(1,M,N,a),dp))/atomic_positions(N,4))
          end do
       end do
 !!      do a=1,3
  !!        mode_eff_charge(a,M) = mode_osc_1(a,M)/sqrt(norm(M))
  !!     end do
       do a=1,3
          do j=1,3
             !         mode_osc(a,j,M) = mode_osc_1(a,M) * mode_osc_1(j,M)
             mode_osc(a,j,M) = mode_osc_1(a,M) * mode_osc_1(j,M) / (one_debye_to_eV**2.0_dp)
             !           write(*,*) a, mode_osc(a,a,M)
          end do
       end do
!       write(*,*) M, phonon_eigenvalues(1,M), mode_osc(1,1,M)+mode_osc(1,2,M)+mode_osc(1,3,M)+mode_osc(2,1,M)+mode_osc(2,2,M)+mode_osc(2,3,M)+mode_osc(3,1,M)+mode_osc(3,2,M)+mode_osc(3,3,M)
       write(*,'(i5,1x,f15.10,6f15.10)') M, phonon_eigenvalues(1,M), mode_osc(1,1,M), mode_osc(2,2,M), mode_osc(3,3,M), mode_osc(2,3,M), mode_osc(3,1,M), mode_osc(1,2,M)
    end do
    
  endsubroutine phonon_eels_calculate_mode_osc
  
  subroutine phonon_eels_calculate_lf_dielectric_tensor(lf_dielectric_tensor, inf_dielectric_tensor, broadening)
    use od_constants, only: pi
    use od_cell, only: cell_volume
    implicit none 

    integer :: nomega
   
    real(dp), intent(inout) :: lf_dielectric_tensor(:, :, :) ! 1:3, 1:3, 1:nomega
    real(dp), intent(in) :: inf_dielectric_tensor(1:3, 1:3) ! epsilon^infinity
    real(dp), intent(in) :: broadening
    complex(dp) :: ctemp
!    real(dp), intent(in) :: born_effective(:,:)

    integer :: gamma_qpoint ! the index for the qpoint at gamma
    integer :: iomega, ir, jr, imode

      
    ! still need
    ! born_effective
    ! iqpoint <-- just a the gamma point

!    gamma_qpoint=phonon_eels_find_gamma_qpoint ! function

!    do iomega=1,nomega
!           lf_dielectric_tensor(ir,jr,iomega)=inf_dielectric_tensor(ir,jr)
!    end do 
           
!    do iomega=1,nomega
!      do ir=1,3
!        do jr=1,3
!          do imode=1,num_eigenvalues
!            ! is this the right way to do i*broadening
!            ctemp = born_effective(imode,ir)*born_effective(imode,jr)/(phonon_eigenvalues(gamma_qpoint,imode)**2- (freq_scale(iomega)+(0, 1)*broadening)**2)
!          end do
!          lf_dielectric_tensor(ir,jr,iomega)=inf_dielectric_tensor(ir,jr)+ 4*pi*ctemp/cell_volume
!        end do
!      end do
!    end do

  end subroutine phonon_eels_calculate_lf_dielectric_tensor

  complex function polarizability(ibin, kx, ky)
    implicit none
    real, intent(in) ::  kx, ky
    integer, intent(in):: ibin

    real(dp) :: K ! Normalisation
    real(dp) :: rtemp ! temporary real

    K=sqrt(kx**2+ky**2)

 !   rtemp = sqrt(lf_dielectric_tensor(3:3,ibin)(kx**2*lf_dielectric_tensor(1:1,ibin) + ky**2*lf_dielectric_tensor(2:2,ibin))/ K**2)

    polarizability= (rtemp-1.0_dp)/(rtemp+1.0_dp)

  end function polarizability

   subroutine phonon_eels_read_phonon_file
    use od_io, only: stdout, maxlen, io_file_unit, io_error
    use od_io, only: seedname
    use od_parameters, only: iprint
    implicit none

    integer :: phonon_in_unit, ierr, i, iatom, iqpoint, ieigenvalue
    integer :: iqpoint_dummy, idummy
    character(maxlen) :: dummy

    character(len=maxlen), allocatable :: atom_name(:)
    real(dp), allocatable :: qpoint_weights(:)
    real(dp) :: rdummy(1:6)

    if (iprint > 2) write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'

    if (allocated(phonon_eigenvectors)) return

    phonon_in_unit = io_file_unit()

    open (unit=phonon_in_unit, file=trim(seedname)//".phonon", form='formatted', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open .phonon file in phonon_eels_read_phonon_file")

    read (phonon_in_unit, *) dummy ! BEGIN header
    read (phonon_in_unit, *) dummy, dummy, dummy, num_ions
    read (phonon_in_unit, *) dummy, dummy, dummy, num_eigenvalues
    read (phonon_in_unit, *) dummy, dummy, dummy, num_qpoints
    read (phonon_in_unit, '(a)') dummy ! Frequency Units
    read (phonon_in_unit, '(a)') dummy ! IR units
    read (phonon_in_unit, '(a)') dummy ! Raman Units
    read (phonon_in_unit, '(a)') dummy ! Unit Cell Vs
    read (phonon_in_unit, *) (phonon_lattice(1, i), i=1, 3) ! Is this the correct way round?
    read (phonon_in_unit, *) (phonon_lattice(2, i), i=1, 3)
    read (phonon_in_unit, *) (phonon_lattice(3, i), i=1, 3)
    read (phonon_in_unit, '(a)') dummy ! Fractional Coodinates

    if (iprint > 2) then
      write (stdout, '(1x,a1,a26,1x,i5,10x,a1)') "|", " Number of Ions: ", num_ions, "|"
      write (stdout, '(1x,a1,a26,1x,i5,10x,a1)') "|", " Number of eigenvalues: ", num_eigenvalues, "|"
      write (stdout, '(1x,a1,a26,1x,i5,10x,a1)') "|", " Number of q-points:", num_qpoints, "|"
    endif

    allocate (atomic_positions(num_ions, 1:4), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate atomic_positions")
    allocate (atom_name(num_ions), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate atom_name")
    allocate (phonon_eigenvectors(num_qpoints, num_eigenvalues, num_ions, 1:3), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate phonon_eigenvectors")
    allocate (phonon_eigenvalues(num_qpoints, num_eigenvalues), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate num_qpoints")
    allocate (qpoint_positions(num_qpoints, 3), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate qpoint_positions")
    allocate (qpoint_weights(num_qpoints), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate qpoint_weights")
    atomic_positions = 0.0_dp
    
    do iatom = 1, num_ions
      read (phonon_in_unit, *) dummy, (atomic_positions(iatom, i), i=1, 3), atom_name(iatom), atomic_positions(iatom, 4)
    enddo
    read (phonon_in_unit, '(a)') dummy ! END header

    do iqpoint = 1, num_qpoints ! Loop over k
       read (phonon_in_unit, *) dummy, iqpoint_dummy, (qpoint_positions(iqpoint, i), i=1, 3), qpoint_weights(iqpoint)
 !      write(*,*) iqpoint_dummy, iqpoint
      if (iqpoint_dummy .ne. iqpoint) call io_error(" ERROR: Error reading q-pt in phonon_eels_read_phonon_file") ! LOTO splitting messes this up 
      do ieigenvalue = 1, num_eigenvalues
        read (phonon_in_unit, *) dummy, phonon_eigenvalues(iqpoint, ieigenvalue)
      enddo
      read (phonon_in_unit, '(a)') dummy ! Phonon Eigenvectors
      read (phonon_in_unit, '(a)') dummy ! Mode, Ion, X, Y, Z

      ! loop over q
      do ieigenvalue = 1, num_eigenvalues
        do iatom = 1, num_ions
          ! Want to do this, but the file is not written in Fortran formatted complex numbers
          !read(phonon_in_unit,*) idummy, idummy, (phonon_eigenvectors(iqpoint, ieigenvalue, iatom, i), i=1,3)
          read (phonon_in_unit, *) idummy, idummy, rdummy(1:6)
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 1) = cmplx(rdummy(1), rdummy(2), dp)
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 2) = cmplx(rdummy(3), rdummy(4), dp)
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 3) = cmplx(rdummy(5), rdummy(6), dp)
        enddo
      enddo
    enddo

    close (phonon_in_unit)

  end subroutine phonon_eels_read_phonon_file

 
  subroutine phonon_eels_read_task
    use od_io, only: stdout, io_error
    use od_parameters, only: phonon_eels_task
    implicit none

    selectcase (phonon_eels_task)
    case ("aloof")
      aloof_scattering = .true.
    case ("impact")
      impact_scattering = .true.
    case ("all")
      impact_scattering = .true.
      aloof_scattering = .true.
    case ('')
      call io_error(' phonon_eels_read_task : No phonon_eels_task found')
    case default
      call io_error(' phonon_eels_calculate : Cannot read phonon_eels_task.')
    end select

  end subroutine phonon_eels_read_task

  subroutine phonon_eels_read_aloof_method
    use od_io, only: stdout, io_error
    use od_parameters, only: phonon_eels_aloof_method
    ! use od_parameters, only: phonon_eels_aloof_methof
    implicit none

    write (stdout, '(a30,a30)') "phonon_eels_aloof_method: ",phonon_eels_aloof_method

    selectcase (phonon_eels_aloof_method)
    case ("dipole")
      dipole_aloof = .true.
    case ("semiclassical")
      semiclassical_aloof = .true.
    case ("all")
      dipole_aloof = .true.
      semiclassical_aloof = .true.
    case ('')
      call io_error(' phonon_eels_read_task : No phonon_eels_aloof_method found')
    case default
      call io_error(' phonon_eels_calculate : Cannot phonon_eels_aloof_method.')
    end select

  end subroutine phonon_eels_read_aloof_method

end module od_phonon_eels
