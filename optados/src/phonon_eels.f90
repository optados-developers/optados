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

  real(kind=dp), allocatable, save  :: adf(:, :)
  real(kind=dp), allocatable, save  :: debye_waller(:)
  complex(kind=dp), allocatable, save  :: phonon_eigenvectors(:, :, :, :) ! iqpoint, ieigenvalues, iatom, i=1,3)
  real(kind=dp), allocatable, save  :: phonon_eigenvalues(:, :)

  real(kind=dp), save  :: phonon_lattice(1:3, 1:3) ! Maybe this is worth keeping
  integer :: num_ions ! found in phonon file
  integer :: num_eigenvalues
  integer :: num_qpoints

contains

  subroutine phonon_eels_calculate
    use od_io, only: stdout, io_error
    use od_electronic, only: elec_read_optical_mat, elec_read_band_energy
    use od_parameters, only: iprint, phonon_eels_task, phonon_eels_aloof_method

    implicit none

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
    call phonon_eels_read_task

    if (impact_scattering) then
      !  call phonon_eels_get_thermal_noise()
      !    call elec_read_band_energy()
      call phonon_eels_read_phonon_file()
      !    call elec_read_optical_mat()
      call phonon_eels_read_chge_trans()
    endif

    if (aloof_scattering) then
      call phonon_eels_read_aloof_method
      if (dipole_aloof) then
        call io_error(' phonon_eels_calculate : Dipole Aloof Not yet implememnted.')
      elseif (semiclassical_aloof) then
        call phonon_eels_read_phonon_file()
        call phonon_eels_read_chge_trans()
      else
        call io_error(' phonon_eels_calculate : Unknown aloof method.')
      endif
    endif

    ! Do the Maths

    ! Write it out

  end subroutine phonon_eels_calculate

  subroutine phonon_eels_get_thermal_noise
    use od_io, only: stdout, io_error
    use od_parameters, only: thermal_noise
    use od_constants, only: preset_debye_waller
    implicit none
    write (stdout, *) "Read phonon_eels_get_thermal_noise()"

    selectcase (thermal_noise)
    case ("adf")
      call phonon_eels_read_adf()
      call phonon_eels_make_debwall()
    case ("debwall", "debye-waller")
      call phonon_eels_read_debwall()
    case ("internal", "optados")
      ! loop over atoms do we know that from phonon?
      ! if atom type == preset_debye_waller then
      ! debeye_waller=preset_debye_waller
      ! endif
    case ('')
      call io_error(' phonon_eels_get_thermal_noise : No thermal_noise found.')
    case default
      call io_error(' phonon_eels_get_thermal_noise : Cannot read thermal_noise.')
    end select
  endsubroutine phonon_eels_get_thermal_noise

  subroutine phonon_eels_read_adf
    use od_constants, only: dp
    use od_parameters, only: iprint
    use od_io, only: stdout, maxlen, seedname, io_file_unit, io_error
    implicit none

    integer :: adf_in_unit, nlines, iline, ierr, i, natoms
    character(len=maxlen) :: dummy, dummy2
    character(len=maxlen), allocatable :: atom_name(:)

    real(dp) :: temperature

    adf_in_unit = io_file_unit()

    open (unit=adf_in_unit, file=trim(seedname)//".adf", form='formatted', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open .adf file in phonon_eels_read_adf")

    ! Go through the file once, check the number of lines
    nlines = 0
    do
      read (adf_in_unit, *, iostat=ierr)
      if (ierr /= 0) exit
      nlines = nlines + 1
    end do
    rewind (unit=adf_in_unit)

    ! Allocate the arrays associated with the file
    if (.not. allocated(adf)) then
      allocate (adf(1:nlines - 4, 1:6), stat=ierr)
    endif
    if (ierr /= 0) call io_error(" Error : cannot allocate adf")
    if (.not. allocated(atom_name)) then
      allocate (atom_name(nlines - 4), stat=ierr)
    endif
    if (ierr /= 0) call io_error(" Error : cannot allocate atom_name")

    ! Read off the comments
    if (iprint > 2) write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    do
      read (adf_in_unit, '(A)') dummy
      if (index(dummy, '#') == 0) then
        ! we do not have a comment line
        exit
      end if
      if (iprint > 2) write (stdout, *) dummy
    end do

    read (dummy, *) dummy2, dummy2, dummy2, natoms

    read (adf_in_unit, *) dummy, temperature

    if (iprint > 2) then
      write (stdout, '(1x,a1,a6,1x,i5,30x,a1)') "|", " natoms:", natoms, "|"
      write (stdout, '(1x,a1,a15,1x,f10.2,30x,a1)') "|", " temperature:", temperature, "|"
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    endif

    if (natoms .ne. nlines - 4) call io_error("Cannot find the number of atoms reported in .adf file.")

    do iline = 1, nlines - 4
      ! would be better if we read in atom_name(iline) number_in_species(iline)
      read (adf_in_unit, *) atom_name(iline), (adf(iline, i), i=1, 6)
    end do

    close (adf_in_unit)

    if (iprint > 2) then
      write (stdout, '(1x,a1,24x,a27,25x,a1)') "|", "Atomic Displacement Factors", "|"
      write (stdout, '(1x,a1,3x,1x,6a11,6x,a1)') "|", "U11", "U22", "U33", "U23", "U31", "U12", "|"
      do iline = 1, nlines - 4
        write (stdout, '(1x,a1,3x,a5,1x,6(e10.3,1x),1x,a1)') "|", trim(atom_name(iline)), (adf(iline, i), i=1, 6), "|"
      enddo
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    endif

  endsubroutine phonon_eels_read_adf

  subroutine phonon_eels_make_debwall
    use od_io, only: stdout
    implicit none
    write (stdout, *) "Read phonon_eels_make_debwall()"
  endsubroutine phonon_eels_make_debwall

  subroutine phonon_eels_read_debwall
    use od_io, only: stdout
    implicit none
    write (stdout, *) "Read phonon_eels_read_debwall()"
  endsubroutine phonon_eels_read_debwall

  subroutine phonon_eels_read_phonon_file
    use od_io, only: stdout, maxlen, io_file_unit, io_error
    use od_io, only: seedname
    use od_parameters, only: iprint
    implicit none

    integer :: phonon_in_unit, ierr, i, iatom, iqpoint, ieigenvalue
    integer :: iqpoint_dummy, idummy
    character(maxlen) :: dummy

    real(dp), allocatable :: atomic_positions(:, :) ! one day this ought to go in the right module.
    character(len=maxlen), allocatable :: atom_name(:)
    real(dp), allocatable :: qpoint_positions(:, :)
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

    allocate (atomic_positions(num_ions, 1:3), stat=ierr)
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

    do iatom = 1, num_ions
      read (phonon_in_unit, *) dummy, (atomic_positions(iatom, i), i=1, 3), atom_name(iatom), dummy
    enddo
    read (phonon_in_unit, '(a)') dummy ! END header

    do iqpoint = 1, num_qpoints ! Loop over k
      read (phonon_in_unit, *) dummy, iqpoint_dummy, (qpoint_positions(iqpoint, i), i=1, 3), qpoint_weights(iqpoint)
      if (iqpoint_dummy .ne. iqpoint) call io_error(" ERROR: Error reading q-pt in phonon_eels_read_phonon_file")
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
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 1) = complex(rdummy(1), rdummy(2))
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 2) = complex(rdummy(3), rdummy(4))
          phonon_eigenvectors(iqpoint, ieigenvalue, iatom, 3) = complex(rdummy(5), rdummy(6))
        enddo
      enddo
    enddo

    close (phonon_in_unit)

  end subroutine phonon_eels_read_phonon_file

  subroutine phonon_eels_read_chge_trans
    use od_io, only: stdout
    implicit none
    write (stdout, *) "Read phonon_eels_read_chge_trans()"
  endsubroutine phonon_eels_read_chge_trans

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

    write (stdout, *) "phonon_eels_aloof_method: "!, phonon_eels_aloof_method

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
