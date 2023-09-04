!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!
! This module contains GPL routines from Wannier90           !
! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          !
!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This version (c) Jonathan Yates 2010                       !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `COPYING' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!
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
module od_io

  use od_constants, only: dp

  implicit none

  private

  integer, public, save           :: stdout
  integer, public, save           :: stderr
  integer, parameter, public :: filename_len = 80
  character(len=filename_len), public, save :: seedname
  character(len=filename_len), public, save :: options
  character(len=filename_len), public, save :: temp_dir
  integer, parameter, public :: maxlen = 120  ! Max column width of input file

  public :: io_get_seedname
  public :: io_time
  public :: io_date
  public :: io_error
  public :: io_file_unit

contains

  subroutine io_get_seedname()
    !==================================================================!
    !                                                                  !
    ! Get the seedname from the commandline                            !
    ! Note iargc and getarg are not standard                           !
    ! Some platforms require them to be external or provide            !
    ! equivalent routines. Not a problem in f2003!                     !
    !===================================================================

    implicit none

    integer :: num_arg

    num_arg = command_argument_count()
    if (num_arg == 0) then
      seedname = '--help'
    elseif (num_arg == 1) then
      call get_command_argument(1, seedname)
      ! Added by F. Mildner to allow for multi_output runs
    elseif (num_arg == 3) then
      call get_command_argument(1, options)
      call get_command_argument(2, temp_dir)
      call get_command_argument(3, seedname)
    else
      call get_command_argument(1, seedname)
      !do something else
    end if

    ! If on the command line the whole seedname.odi was passed, I strip the last ".win"
    if (len(trim(seedname)) .ge. 5) then
      if (seedname(len(trim(seedname)) - 4 + 1:) .eq. ".odi") then
        seedname = seedname(:len(trim(seedname)) - 4)
      end if
    end if

  end subroutine io_get_seedname

  !==================================================================!
  subroutine io_error(error_msg)
    !==================================================================!
    !                                                                  !
    ! Aborts giving error message                                      !
    !                                                                  !
    !===================================================================

    implicit none
    character(len=*), intent(in) :: error_msg
    write (stderr, *) 'Exiting.......'
    write (stderr, '(1x,a)') trim(error_msg)
    close (stderr)

    error stop "Optados error: examine the output/error file for details"

  end subroutine io_error

  !==================================================================!
  subroutine io_date(cdate, ctime)
    !==================================================================!
    !                                                                  !
    !     Returns two strings containing the date and the time         !
    !     in human-readable format. Uses a standard f90 call.          !
    !                                                                  !
    !===================================================================
    implicit none
    character(len=11), intent(out) :: cdate
    character(len=9), intent(out) :: ctime

    character(len=3), dimension(12) :: months
    data months/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
    integer date_time(8)
    !
    call date_and_time(values=date_time)
    !
    write (cdate, '(i2,1x,a3,1x,i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime, '(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)

  end subroutine io_date

  !==================================================================!
  function io_time()
    !==================================================================!
    !                                                                  !
    ! Returns elapsed CPU time in seconds since its first call         !
    ! uses standard f90 call                                           !
    !                                                                  !
    !===================================================================
    use od_constants, only: dp
    implicit none

    real(kind=dp) :: io_time

    ! t0 contains the time of the first call
    ! t1 contains the present time
    real(kind=dp) :: t0, t1
    logical :: first = .true.
    save first, t0
    !
    call cpu_time(t1)
    !
    if (first) then
      t0 = t1
      io_time = 0.0_dp
      first = .false.
    else
      io_time = t1 - t0
    end if
    return
  end function io_time

  !==================================================================!
  function io_file_unit()
    !==================================================================!
    !                                                                  !
    ! Returns an unused unit number                                    !
    ! (so we can open a file on that unit                              !
    !                                                                  !
    !===================================================================
    implicit none

    integer :: io_file_unit, unit
    logical :: file_open

    unit = 9
    file_open = .true.
    do while (file_open)
      unit = unit + 1
      inquire (unit, OPENED=file_open)
    end do

    io_file_unit = unit

    return
  end function io_file_unit

end module od_io
