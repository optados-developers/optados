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
program optados
  !=========================================================================!
  !                             O P T A D O S                               !
  !                      OPTics And Density Of States                       !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  ! Described below.                                                        !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  !-------------------------------------------------------------------------!
  ! Written by Andrew Morris, Rebecca Nicholls, Chris Pickard               !
  !             and Jonathan Yates      2010                                !
  !=========================================================================!
  use od_comms, only: comms_setup, on_root, comms_end, num_nodes
  use od_constants, only: dp
  use od_io, only: io_get_seedname, io_time, io_date, io_file_unit,&! Functions
       & stdout, stderr, seedname, options, temp_dir                ! Variables
  use od_parameters, only: param_read, param_write_header, param_Dist, param_write, &
    param_dealloc, pdos, pdis, dos, jdos, core, optics, photo, iprint, param_write_atomic_coord, &
    devel_flag, photo_photon_energy, photo_model
  use od_cell, only: cell_calc_lattice, cell_report_parameters, cell_dist
  use od_electronic, only: elec_read_band_energy, elec_read_band_energy_ordered, elec_report_parameters
  use od_dos, only: dos_calculate
  use od_jdos, only: jdos_calculate
  use od_core, only: core_calculate
  use od_pdos, only: pdos_calculate
  use od_pdis, only: pdis_calculate
  use od_optics, only: optics_calculate
  use od_photo, only: photo_calculate
  use od_build, only: build_info
  implicit none

  real(kind=dp)    :: time0, time1       ! Variables for timing
  logical          :: odo_found         ! Ouptut file exists?
  character(len=9) :: stat, pos          ! Status and position of .odo file
  character(len=9) :: ctime             ! Temp. time string
  character(len=11):: cdate             ! Temp. date string
  character(len=120):: filename            ! Added by Felix Mildner, 03/23 for multi file output

  ! call sleep(25)

  time0 = io_time()

  call comms_setup

  if (on_root) then
    call io_get_seedname()
    ! If blank set to seedname='--help'
    if (trim(seedname) == '-h' .or. trim(seedname) == '--help') call help_output
    if (trim(seedname) == '-v' .or. trim(seedname) == '--version') call version_output
    !-------------------------------------------------------------------------!
    ! O R G A N I S E   T H E   E R R O R   F I L E
    stderr = io_file_unit()
    ! This is to allow multiple OptaDOS photoemission runs to be performed in the same directory.
    if ((index(options, '-temp') > 0)) then
      filename = trim(adjustl(temp_dir))//'/'//trim(seedname)//'.opt_err'
      open (unit=stderr, file=filename)
    else
      open (unit=stderr, file=trim(seedname)//'.opt_err')
    end if
    call io_date(cdate, ctime)
    write (stderr, *) 'OptaDOS: Execution started on ', cdate, ' at ', ctime
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    ! O R G A N I S E   T H E   O U T P U T   F I L E  A N D
    ! R E A D   A N D   W R I T E   U S E R   P A R A M E T E R S
    call param_read()
    ! This is to allow multiple simultaneous OptaDOS photoemission runs to be performed in the same directory.
    if ((index(options, '-temp') /= 0)) then
      filename = trim(adjustl(temp_dir))//'/'//trim(seedname)//'.odo'
      inquire (file=filename, exist=odo_found)
    else
      inquire (file=trim(seedname)//'.odo', exist=odo_found)
    end if
    if (odo_found) then
      stat = 'old'
    else
      stat = 'replace'
    end if
    pos = 'append'

    stdout = io_file_unit()
    ! This is to allow multiple OptaDOS photoemission runs to be performed in the same directory.
    if ((index(options, '-temp') /= 0)) then
      open (unit=stdout, file=filename, status=trim(stat), position=trim(pos))
    else
      open (unit=stdout, file=trim(seedname)//'.odo', status=trim(stat), position=trim(pos))
    end if
    write (stdout, *) 'OptaDOS: Execution started on ', cdate, ' at ', ctime
    write (stdout, '(1x,a26,i5,a10)') 'Parallelised over', num_nodes, ' thread(s)'
    if (iprint > 0) call param_write_header()
    if (iprint > 0) call param_write()
    time1 = io_time()

    if (iprint > 1) write (stdout, '(1x,a40,f11.3,a)') 'Time to read parameters ', time1 - time0, ' (sec)'
    !-------------------------------------------------------------------------!
  end if

  if (pdis) then
    call elec_read_band_energy_ordered
  else
    call elec_read_band_energy
  end if

  if (on_root) then
    call cell_calc_lattice
    if (iprint > 0) call param_write_atomic_coord
    if (iprint > 0) call cell_report_parameters
    if (iprint > 0) call elec_report_parameters
  end if
  ! now send the data from the parameter file to each node

  call param_dist
  call cell_dist

  !-------------------------------------------------------------------------!
  ! C A L L   P D O S   R O U T I N E S
  if (pdos) then
    time0 = io_time()
    call pdos_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to calculate Projected Density of States          &
           &      ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   P D I S   R O U T I N E S
  if (pdis) then
    time0 = io_time()
    call pdis_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to calculate Projected Dispersion &
           &      ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   C O R E   R O U T I N E S
  if (core) then
    time0 = io_time()
    call core_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '|                                                                            |'
      write (stdout, '(1x,a59,f11.3,a8)') &
                       '+ Time to calculate Core Level Spectra       &
                       &                  ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   D O S  R O U T I N E S
  if (dos) then
    time0 = io_time()
    call dos_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a59,f11.3,a8)') &
                       '+ Time to calculate Density of States        &
                       &                  ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   O P T I C S   R O U T I N E S
  if (optics) then
    time0 = io_time()
    call optics_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '|                                                                            |'
      write (stdout, '(1x,a59,f11.3,a8)') &
        '+ Time to calculate Optical properties                         ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   P H O T O E M I S S I O N   R O U T I N E S
  if (photo) then
    time0 = io_time()
    call photo_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '|                                                                            |'
      write (stdout, '(1x,a59,f11.3,a8)') &
        '+ Time to calculate Photoemission                              ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if
  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! C A L L   J D O S   R O U T I N E S
  if (jdos) then
    time0 = io_time()
    call jdos_calculate
    time1 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '|                                                                            |'
      write (stdout, '(1x,a59,f11.3,a8)') &
        '+ Time to calculate Joint Density of States                    ', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, *) ' '
      write (stdout, *) ' '
    end if

  end if
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! F I N A L I S E
  call param_dealloc

  if (on_root) then
    call io_date(cdate, ctime)
    write (stdout, *)
    write (stdout, *) 'OptaDOS: Execution complete on ', cdate, ' at ', ctime

    close (stdout)
    close (stderr, status='delete')
  end if

  call comms_end

  !-------------------------------------------------------------------------!
contains
  subroutine help_output
    use od_constants, only: optados_version, copyright
    implicit none
    write (*, *)
    write (*, *) " OptaDOS version ", trim(build_info%build)
    write (*, *)
    write (*, *) " Andrew J. Morris, R. J. Nicholls, C. J. Pickard and J. R. Yates", trim(copyright)
    write (*, *) " Usage: optados <seedname>"
    write (*, *)
    stop
  end subroutine help_output

  subroutine version_output
    use od_build, only: build_info
    use od_constants, only: optados_version, copyright
    implicit none
    write (*, *)
    write (*, *) " OptaDOS version ", trim(build_info%build)
    write (*, *)
    write (*, *) " Andrew J. Morris, R. J. Nicholls, C. J. Pickard and J. R. Yates", trim(copyright)
    write (*, *) " Compiled with "//trim(build_info%compiler)//" on "//trim(build_info%compile_date)&
         & //" at "//trim(build_info%compile_time)//"."
    write (*, *) " Compile type: "//trim(build_info%build_type)//", "//trim(build_info%comms_arch)
    write (*, *) " From source "//trim(build_info%build)//" submitted on "//trim(build_info%source_date)&
         &//" at "//trim(build_info%source_time)//"."

    stop
  end subroutine version_output

end program optados
