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
!===========================================================================!
! MODULE od_read_odd                                                        !
!                                                                           !
! J. A. J. Whaley-Baldwin, May 2026                                         !
!                                                                           !
! This module contains subroutines that parse data from an OptaDOS .odd     !
! file. Functionality is provided to parse-in individual blocks, as well    !
! as wrapper subroutines that parse all quantities relevant for a certain   !
! task (e.g. for phonon_eels).                                              !
!===========================================================================!
module od_read_odd

  ! OptaDOS IO.
  use od_io, only: seedname, stdout, io_file_unit, io_error

  ! OptaDOS constants.
  use od_constants, only: dp

  ! OptaDOS parameters.
  use od_parameters, only: iprint, phonon_eels_task

  implicit none

  private

  ! Exposed routines.
  public :: odd_read_vib_eels_data
  public :: odd_read_inf_dielectric_tensor
  public :: odd_read_born_eff_charges
  public :: odd_read_partial_charges
  public :: odd_read_atomic_displacement_params

  ! These are parsed in from the .odd file.
  real(kind=dp),save,public                 :: ADP_temperature                         ! Temperature of supplied ADPs
  real(kind=dp),allocatable,save,public     :: inf_dielectric_tensor(:,:)              ! (i,j)
  real(kind=dp),allocatable,save,public     :: born_eff_ch_tensor(:,:,:)               ! (atom_idx,i,j)
  real(kind=dp),allocatable,save,public     :: partial_charges(:)                      ! (atom_idx)
  real(kind=dp),allocatable,save,public     :: atomic_displacement_params(:,:,:)       ! (atom_idx,i,j)

  ! Flags to check whether quantities relevant for certain tasks have been parsed or not.
  logical,save,public :: parsed_eels_data_from_odd = .false.

  contains

    !=========================================================================!
    subroutine odd_read_vib_eels_data
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                       !
    !                                                                         !
    ! This subroutine is a wrapper, that simply calls all of the parsing      !
    ! subroutines in this module that are relevant for the requested          !
    ! phonon EELS task.                                                       !
    !                                                                         !
    ! This should be called before any vib-EELS subroutines are used (if      !
    ! not, the vib-EELS subroutines will call it anyway).                     !
    !=========================================================================!
      implicit none

      if ( index(phonon_eels_task, 'impact') > 0 ) then
        CALL odd_read_partial_charges
        CALL odd_read_atomic_displacement_params
      else if ( index(phonon_eels_task, 'aloof') > 0 ) then
        CALL odd_read_inf_dielectric_tensor
        CALL odd_read_born_eff_charges
      else if ( index(phonon_eels_task, 'all') > 0 ) then
        CALL odd_read_partial_charges
        CALL odd_read_atomic_displacement_params
        CALL odd_read_inf_dielectric_tensor
        CALL odd_read_born_eff_charges

      else
        CALL io_error("ERROR: phonon_eels_task not recognized. Should be one of 'impact', 'aloof', or 'all'")

      end if

      parsed_eels_data_from_odd = .true.

    end subroutine odd_read_vib_eels_data

    !=========================================================================!
    subroutine odd_read_inf_dielectric_tensor
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                     !
    !                                                                         !
    ! This subroutine reads in a .odd file, and parses the                    !
    ! (zero frequency) electronic dielectric tensor (EPSILON_INF block).      !
    !                                                                         !
    ! Calling this will allocate & populate the 'inf_dielectric_tensor'       !
    ! array.                                                                  !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      integer            :: odd_data_unit, ierr, i, parse_stat, row
      character(len=256) :: line
      logical            :: read_epsilon_block
      real(kind=dp)      :: epsilon_inf(3,3)

      ! Allocate eps_inf array.
      allocate(inf_dielectric_tensor(3,3))

      odd_data_unit = io_file_unit()

      read_epsilon_block = .false.
      row = 0

      open(unit=odd_data_unit, file=trim(seedname)//".odd", status="old", action="read")

      do
        read(odd_data_unit,'(A)',iostat=ierr) line
        if (ierr /= 0) exit

        ! Skip any comment lines.
        if (index(adjustl(line), '#') == 1) cycle

        ! Detect start of epsilon_inf block.
        if (trim(adjustl(line)) == "BEGIN EPSILON_INF") then
          read_epsilon_block = .true.
          row = 0
          cycle
        end if

        ! Detect end of epsilon_inf block.
        if (trim(adjustl(line)) == "END EPSILON_INF") then
          read_epsilon_block = .false.
          exit
        end if

        ! Read data inside epsilon_inf block, and handle error cases.
        if (read_epsilon_block) then
          if (row .eq. 3) then
            CALL io_error("ERROR: Too many rows in EPSILON_INF block in .odd file")
          end if
          row = row + 1
          read(line, *, iostat=parse_stat) inf_dielectric_tensor(row, 1:3)
          if (parse_stat /= 0) then
            CALL io_error("ERROR: Problem parsing EPSILON_INF block in .odd file")
          end if
        end if
      end do

      close(odd_data_unit)

      ! For debug: Print inf. dielectric tensor (based on iprint level).
      if (iprint > 2) then
        write(stdout,*) "INF. DIELECTRIC TENSOR (parsed from "//TRIM(seedname)//".odd):"
        write(stdout,*) ""
        do i=1,3
          write(stdout,"(F8.5,F8.5,F8.5)") inf_dielectric_tensor(i,1), inf_dielectric_tensor(i,2), inf_dielectric_tensor(i,3)
        end do
        write(stdout,*) ""
      end if

    end subroutine odd_read_inf_dielectric_tensor

    !=========================================================================!
    subroutine odd_read_born_eff_charges
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                     !
    !                                                                         !
    ! This subroutine reads in a .odd file, and parses the Born effective     !
    ! charge tensors for each atom (BORN_EFF_CHARGE_TENSOR block).            !
    !                                                                         !
    ! Calling this will allocate & populate the 'born_eff_ch_tensor' array.   !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      integer                    :: odd_data_unit, ierr, i, N, atom_count, idx
      character(len=256)         :: line
      logical                    :: read_eff_charge_block
      real(kind=dp)              :: row_values(3)

      odd_data_unit = io_file_unit()

      read_eff_charge_block = .false.
      atom_count = 0

      open(unit=odd_data_unit, file=trim(seedname)//".odd", status="old", action="read")

      ! First pass; count atoms.
      read_eff_charge_block = .false.
      atom_count = 0

      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit

        line = adjustl(line)

        if (len_trim(line) == 0) cycle
        if (index(line, '#') == 1) cycle

        if (trim(line) == "BEGIN BORN_EFF_CHARGE_TENSOR") then
          read_eff_charge_block = .true.
          cycle
        end if

        if (trim(line) == "END BORN_EFF_CHARGE_TENSOR") then
          exit
        end if

        if (read_eff_charge_block) then
          ! Each atom starts with an integer index line.
          read(line, *, iostat=ierr) idx
          if (ierr == 0) atom_count = atom_count + 1
        end if
      end do

      if (atom_count == 0) then
        CALL io_error("ERROR: No Born effective charge tensors found in .odd file")
      end if

      ! Allocate born_eff_ch_tensor.
      allocate(born_eff_ch_tensor(atom_count,3,3))

      ! Second pass; fill data.
      rewind(odd_data_unit)
      read_eff_charge_block = .false.

      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit

        line = adjustl(line)

        if (len_trim(line) == 0) cycle
        if (index(line, '#') == 1) cycle

        if (trim(line) == "BEGIN BORN_EFF_CHARGE_TENSOR") then
          read_eff_charge_block = .true.
          atom_count = 0
          cycle
        end if

        if (trim(line) == "END BORN_EFF_CHARGE_TENSOR") then
          exit
        end if

        if (read_eff_charge_block) then

          ! Read atom index.
          read(line, *, iostat=ierr) idx
          if (ierr /= 0) then
            CALL io_error("ERROR: Problem parsing BORN_EFF_CHARGE_TENSOR block in .odd file")
          end if

          atom_count = atom_count + 1

          ! Read 3x3 tensor.
          do i = 1,3
            read(odd_data_unit, '(A)', iostat=ierr) line
            if (ierr /= 0) then
              CALL io_error("ERROR: Unexpected EOF inside BEC tensor within .odd file")
            end if

            line = adjustl(line)
            read(line, *, iostat=ierr) row_values

            if (ierr /= 0) then
              CALL io_error("ERROR: Problem parsing BORN_EFF_CHARGE_TENSOR block in .odd file")
            end if

            born_eff_ch_tensor(atom_count,i,1:3) = row_values
          end do

        end if
      end do

      close(odd_data_unit)

      ! For debug: Print BEC tensor (based on iprint level).
      if (iprint > 2) then
        write(stdout,*) "BORN EFFECTIVE CHARGE TENSOR (parsed from "//TRIM(seedname)//".odd):"
        write(stdout,*) ""
        do N=1,atom_count
          write(stdout,*) "Atom ",N
          do i=1,3
            write(stdout,"(F8.5,F8.5,F8.5)") born_eff_ch_tensor(N,i,1), born_eff_ch_tensor(N,i,2), born_eff_ch_tensor(N,i,3)
          end do
          write(stdout,*) ""
        end do
        write(stdout,*) ""
      end if

    end subroutine odd_read_born_eff_charges

    !=========================================================================!
    subroutine odd_read_partial_charges
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                     !
    !                                                                         !
    ! This subroutine reads in a .odd file, and parses the partial charges on !
    ! each atom (PARTIAL_CHARGES block).                                      !
    !                                                                         !
    ! Calling this will allocate & populate the 'partial_charges' array.      !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      character(len=256) :: line
      integer            :: odd_data_unit, ierr
      logical            :: read_partial_charge_block
      integer            :: idx, atom_count
      real(kind=dp)      :: charge, init_check_value

      odd_data_unit = io_file_unit()
      read_partial_charge_block = .false.

      open(unit=odd_data_unit, file=trim(seedname)//".odd", status="old", action="read")

      ! First pass; count number of atoms.
      atom_count = 0
      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        ! Check for start of partial_charge block.
        if (index(line, "BEGIN PARTIAL_CHARGES") > 0) then
          read_partial_charge_block = .true.
          cycle
        end if

        ! Increment atom_count.
        if (read_partial_charge_block) then
          atom_count = atom_count + 1
        end if

        ! Check for end of partial_charge block.
        if (index(line, "END PARTIAL_CHARGES") > 0) then
          read_partial_charge_block = .false.
          exit
        end if
      end do
      atom_count = atom_count - 1

      ! Allocate partial_charges.
      allocate(partial_charges(atom_count))

      ! Initially set all partial_charges to a nonsensical value.
      ! This is used to check that the data has been read-in correctly, after the loop.
      init_check_value = 1.0E6_dp
      partial_charges = init_check_value

      ! Second pass; parse charges and populate partial_charges.
      rewind(odd_data_unit)
      read_partial_charge_block = .false.
      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        ! Check for start of partial_charge block.
        if (index(line, "BEGIN PARTIAL_CHARGES") > 0) then
          read_partial_charge_block = .true.
          cycle
        end if

        ! Check for end of partial_charge block.
        if (index(line, "END PARTIAL_CHARGES") > 0) then
          read_partial_charge_block = .false.
          exit
        end if

        ! If inside the block, parse partial_charge.
        if (read_partial_charge_block) then
          if (len_trim(line) == 0) cycle
          read(line, *, iostat=ierr) idx, charge
          if (ierr == 0) then
          partial_charges(idx) = charge
          else
            CALL io_error("ERROR: Problem parsing PARTIAL_CHARGES block in .odd file")
          end if
        end if
      end do

      close(odd_data_unit)

      ! Check that the partial charges have been read-in correctly.
      ! If so, then the following condition will not be met.
      if ( ANY(ABS(partial_charges - init_check_value) < 1.0E-3) ) then
        CALL io_error("ERROR: One or more partial charges is missing in the .odd file, and/or has been defined multiple times")
      end if

      ! For debug: Print partial charges (based on iprint level).
      if (iprint > 2) then
        write(stdout,*) "PARTIAL CHARGES (parsed from "//TRIM(seedname)//".odd):"
        write(stdout,*) ""
        do idx=1,atom_count
          write(stdout,*) partial_charges(idx)
        end do
        write(stdout,*) ""
      end if

    end subroutine odd_read_partial_charges

    !=========================================================================!
    subroutine odd_read_atomic_displacement_params
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                     !
    !                                                                         !
    ! This subroutine reads in a .odd file, and parses the atomic             !
    ! displacement parameters (ADPs) for each atom, alongside the temperature !
    ! at which the ADPs were supplied.                                        !
    !                                                                         !
    ! Calling this will allocate & populate the 'atomic_displacement_params'  !
    ! array, and will set 'ADP_temperature'.                                  !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      integer            :: odd_data_unit, ierr, i, idx, atom_count
      real(kind=dp)      :: dummy, init_check_value, Uxx, Uyy, Uzz, Uyz, Uzx, Uxy
      character(len=256) :: line
      logical            :: read_atomic_displacement_params_block,expecting_temperature

      odd_data_unit = io_file_unit()
      open(unit=odd_data_unit, file=trim(seedname)//".odd", status="old", action="read")
      read_atomic_displacement_params_block = .false.
      expecting_temperature = .false.

      ! First pass; count atoms.
      atom_count = 0
      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit

        line = adjustl(line)

        if (len_trim(line) == 0) cycle
        if (index(line, '#') == 1) cycle

        if (trim(line) == "BEGIN ATOMIC_DISPLACEMENT_PARAMETERS") then
          read_atomic_displacement_params_block = .true.
          expecting_temperature = .true.
          cycle
        end if

        if (trim(line) == "END ATOMIC_DISPLACEMENT_PARAMETERS") then
          exit
        end if

        if (expecting_temperature) then
          read(line, *, iostat=ierr) ADP_temperature, dummy

          if (ierr == 0) then
            CALL io_error("ERROR: ADP temperature line is missing, or has multiple entries, in .odd file")
          else
            read(line, *, iostat=ierr) ADP_temperature
            if (ierr /= 0) then
              CALL io_error("ERROR: Problem parsing ADP temperature in .odd file")
            end if
          end if

          expecting_temperature = .false.
          cycle
        end if

        if (read_atomic_displacement_params_block) then
          ! Each atom starts with an integer index line.
          read(line, *, iostat=ierr) idx, Uxx, Uyy, Uzz, Uyz, Uzx, Uxy
          if (ierr == 0) atom_count = atom_count + 1
        end if
      end do

      if (atom_count == 0) then
        CALL io_error("ERROR: No Atomic Displacement Parameters found in .odd file")
      end if

      ! Allocate atomic_displacement_params.
      allocate(atomic_displacement_params(atom_count,3,3))

      ! Initially set all ADPs to a nonsensical value.
      ! This is used to check that the data has been read-in correctly, after the loop.
      init_check_value = 1.0E6_dp
      atomic_displacement_params = init_check_value

      ! Second pass; parse temperature, ADPs, and populate atomic_displacement_params.
      rewind(odd_data_unit)
      read_atomic_displacement_params_block = .false.
      expecting_temperature = .false.

      do
        read(odd_data_unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        ! Check for start of ADPs block.
        if (index(line, "BEGIN ATOMIC_DISPLACEMENT_PARAMETERS") > 0) then
          read_atomic_displacement_params_block = .true.
          expecting_temperature = .true.
          cycle
        end if

        ! Check for end of ADPs block.
        if (index(line, "END ATOMIC_DISPLACEMENT_PARAMETERS") > 0) then
          read_atomic_displacement_params_block = .false.
          exit
        end if

        ! If inside the block, parse temperature & ADPs.
        if (read_atomic_displacement_params_block) then
          if (len_trim(line) == 0) cycle
          if (expecting_temperature) then
            read(line, *, iostat=ierr) ADP_temperature, dummy

            if (ierr == 0) then
              CALL io_error("ERROR: Missing ADP temperature line in .odd file")
            else
              read(line, *, iostat=ierr) ADP_temperature
              if (ierr /= 0) then
                CALL io_error("ERROR: Problem parsing ADP temperature in .odd file")
              end if
            end if

            expecting_temperature = .false.
            cycle
          end if
          read(line, *, iostat=ierr) idx, Uxx, Uyy, Uzz, Uyz, Uzx, Uxy
          if (ierr == 0) then
            ! Unique elements.
            atomic_displacement_params(idx,1,1) = Uxx
            atomic_displacement_params(idx,2,2) = Uyy
            atomic_displacement_params(idx,3,3) = Uzz
            atomic_displacement_params(idx,2,3) = Uyz
            atomic_displacement_params(idx,3,1) = Uzx
            atomic_displacement_params(idx,1,2) = Uxy
            ! Duplicates, by symmetry.
            atomic_displacement_params(idx,1,3) = Uzx
            atomic_displacement_params(idx,2,1) = Uxy
            atomic_displacement_params(idx,3,2) = Uyz
          else
            CALL io_error("ERROR: Problem parsing ATOMIC_DISPLACEMENT_PARAMETERS block in .odd file")
          end if
        end if
      end do

      close(odd_data_unit)

      ! Check that the ADPs have been read-in correctly.
      ! If so, then the following condition will not be met.
      if ( ANY(ABS(atomic_displacement_params - init_check_value) < 1.0E-3) ) then
        CALL io_error("ERROR: One or more atomic displacement parameters is missing in the .odd file, &
        and/or has been defined multiple times")
      end if

      ! For debug: Print ADPs (based on iprint level).
      if (iprint > 2) then
        write(stdout,*) "ATOMIC DISPLACEMENT PARAMETERS (parsed from "//TRIM(seedname)//".odd):"
        write(stdout,*) ""
        write(stdout,*) "      U11         U22         U33         U23         U31         U12"
        do idx=1,atom_count
          write(stdout,'(I3,A,F8.6,A,F8.6,A,F8.6,A,F8.6,A,F8.6,A,F8.6)') idx,"    ",atomic_displacement_params(idx,1,1), &
          "    ",atomic_displacement_params(idx,2,2), "    ",atomic_displacement_params(idx,3,3),"    ", &
          atomic_displacement_params(idx,2,3),"    ",atomic_displacement_params(idx,3,1),"    ", &
          atomic_displacement_params(idx,1,2)
        end do
        write(stdout,*) ""
        write(stdout,'(A,F12.6,A)') " TEMPERATURE AT WHICH ADPs WERE SUPPLIED: ",ADP_temperature," K"
        write(stdout,*) ""
      end if

      close(odd_data_unit)

    end subroutine odd_read_atomic_displacement_params

end module od_read_odd