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
      call phonon_eels_get_thermal_noise()
      call elec_read_band_energy()
      call phonon_eels_read_phonon_file()
      call elec_read_optical_mat()
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
    use od_io, only: stdout
    implicit none
    write (stdout, *) "Read phonon_eels_get_thermal_noise()"
  endsubroutine phonon_eels_get_thermal_noise

  subroutine phonon_eels_read_phonon_file
    use od_io, only: stdout
    implicit none
    write (stdout, *) "Read phonon_eels_read_phonon_file()"
  endsubroutine phonon_eels_read_phonon_file

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
    use od_io, only: stdout
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
