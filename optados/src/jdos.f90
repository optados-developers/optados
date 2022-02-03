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
!===============================================================================
module od_jdos
  use od_constants, only: dp
  implicit none

  !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S
  public :: jdos_calculate

contains

  !===============================================================================
  subroutine jdos_calculate
    !===============================================================================
    use od_jdos_utils, only: jdos_utils_calculate, E, jdos_fixed &
         &, jdos_adaptive, jdos_linear
    use od_parameters, only: fixed, adaptive, linear
    use od_io, only: io_time, stdout
    use od_comms, only: on_root
    implicit none

    real(kind=dp) :: time0, time1

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                           Joint Density of States                          +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    call jdos_utils_calculate

    !-------------------------------------------------------------------------------
    ! W R I T E   O U T   J D O S

    time0 = io_time()
    ! Otherwise we have written to wdos and dos, so they can be called
    ! by whatever.
    if (on_root) then
      if (fixed) call write_jdos(E, jdos_fixed, "fixed")
      if (adaptive) call write_jdos(E, jdos_adaptive, "adaptive")
      if (linear) call write_jdos(E, jdos_linear, "linear")
      !if(quad)    call write_jdos(E, dos_quad, intdos_quad, "quad")
    end if
    time1 = io_time()

    !-------------------------------------------------------------------------------

  end subroutine jdos_calculate

  !===============================================================================
  subroutine write_jdos(E, dos, dos_name)
    use od_parameters, only: dos_per_volume, output_format
    use od_jdos_utils, only: jdos_nbins
    use od_electronic, only: nspins
    use od_io, only: io_file_unit, seedname, stdout, io_error, io_date
    !===============================================================================
    ! This routine receives an energy scale, a density of states and a file name
    ! and writes out the DOS to disk
    !===============================================================================
    implicit none
    real(kind=dp), intent(in) :: E(jdos_nbins)
    real(kind=dp), intent(in) :: dos(jdos_nbins, nspins)
    character(len=*), intent(in) :: dos_name
    integer :: i, dos_file, ierr
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=22) :: dos_units, intdos_units

    dos_file = io_file_unit()
    open (unit=dos_file, file=trim(seedname)//'.j'//trim(dos_name)//'.dat', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open output file in dos: write_dos")

    dos_units = "(electrons per eV)"; intdos_units = "(electrons)"
    if (dos_per_volume) then
      dos_units = "(electrons per eV/A^3)"
      intdos_units = "(electrons per A^3)"
    end if

    write (dos_file, *) "##############################################################################"
    write (dos_file, *) "#"
    write (dos_file, *) "#                  O p t a D O S   o u t p u t   f i l e "
    write (dos_file, '(1x,a1)') "#"
    write (dos_file, *) "#    Density of States using ", trim(dos_name), " broadening"
    call io_date(cdate, ctime)
    write (dos_file, *) '#  Generated on ', cdate, ' at ', ctime
    write (dos_file, *) "# Column        Data"
    write (dos_file, *) "#    1        Energy (eV)"
    if (nspins > 1) then
      write (dos_file, *) "#    2        Up-spin DOS ", trim(dos_units)
      write (dos_file, *) "#    3        Down-spin DOS ", trim(dos_units)
      write (dos_file, *) "#    4        Up-spin Integrated DOS ", trim(intdos_units)
      write (dos_file, *) "#    5        Down-spin Integrated DOS ", trim(intdos_units)
    else
      write (dos_file, *) "#    2        DOS ", trim(dos_units)
      write (dos_file, *) "#    3        Integrated DOS ", trim(intdos_units)
    end if
    write (dos_file, '(1x,a1)') "#"
    write (dos_file, '(1x,a78)') "##############################################################################"

    if (nspins > 1) then
      do i = 1, jdos_nbins
        write (dos_file, '(3(E21.13,2x))') E(i), dos(i, 1), -dos(i, 2)
      end do
    else
      do i = 1, jdos_nbins
        write (dos_file, '(2(E21.13,2x))') E(i), dos(i, 1)
      end do
    end if
    close (dos_file)

    if (trim(output_format) == "xmgrace") then
      call write_jdos_xmgrace(dos_name, E, dos)
    elseif (trim(output_format) == "gnuplot") then
      write (stdout, *) " WARNING: GNUPLOT output not yet available, calling xmgrace"
      call write_jdos_xmgrace(dos_name, E, dos)
      !     call write_dos_gnuplot(dos_name,E,dos)
    else
      write (stdout, *) " WARNING: Unknown output format requested, continuing..."
    end if

  end subroutine write_jdos
  !===============================================================================

  !===============================================================================
  subroutine write_jdos_xmgrace(dos_name, E, dos)
    !===============================================================================
    use xmgrace_utils
    use od_jdos_utils, only: jdos_nbins
    use od_electronic, only: nspins
    use od_io, only: io_file_unit, io_error, seedname
    implicit none

    real(kind=dp), intent(in) :: E(jdos_nbins)
    real(kind=dp), intent(in) :: dos(jdos_nbins, nspins)
    real(kind=dp) :: min_x, max_x, min_y, max_y

    integer :: batch_file, ierr
    character(len=*), intent(in) :: dos_name

    batch_file = io_file_unit()
    open (unit=batch_file, file=trim(seedname)//'.j'//trim(dos_name)//'.agr', iostat=ierr)
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open xmgrace batch file in dos: write_jdos_xmgrace")

    min_x = minval(E)
    max_x = maxval(E)

    min_y = 0
    max_y = maxval(dos)
    if (nspins > 1) then
      min_y = -max_y
    end if

    call xmgu_setup(batch_file)
    call xmgu_legend(batch_file)
    call xmgu_title(batch_file, min_x, max_x, min_y, max_y, "Joint Electronic Density of States")
    call xmgu_subtitle(batch_file, "Generated by OptaDOS")

    call xmgu_axis(batch_file, "x", "Energy eV")
    call xmgu_axis(batch_file, "y", "JDOS")

    if (nspins > 1) then
      call xmgu_data_header(batch_file, 0, 1, "up-spin channel")
      call xmgu_data_header(batch_file, 1, 2, "down-spin channel")
      call xmgu_data(batch_file, 0, E(:), dos(:, 1))
      call xmgu_data(batch_file, 1, E(:), -dos(:, 2))
    else
      call xmgu_data_header(batch_file, 0, 1, "Total JDOS")
      call xmgu_data(batch_file, 0, E(:), dos(:, 1))
    end if

    close (batch_file)

  end subroutine write_jdos_xmgrace

end module od_jdos
