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
! E L E C T R O N I C
! Stores variables to do with the electrons and energy eigenvalues in the
! system. At this stage I don't see it contining many functions, but it
! was important for the dos module not to have too many global variables
! AJM Dec 2010
!=========================================================================!

module od_electronic
  use od_constants, only: dp

  implicit none

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  real(kind=dp), allocatable, public, save  :: band_energy(:, :, :)
  real(kind=dp), allocatable, public, save  :: band_gradient(:, :, :, :)
  complex(kind=dp), allocatable, public, save  :: optical_mat(:, :, :, :, :)
  complex(kind=dp), allocatable, public, save  :: elnes_mat(:, :, :, :, :)

  real(kind=dp), public, save :: efermi ! The fermi energy we finally decide on
  logical, public, save       :: efermi_set = .false. ! Have we set efermi?
  real(kind=dp), public, save :: unshifted_efermi ! The fermi energy we finally decide on, perhaps not set to 0
  real(kind=dp), public, save :: efermi_castep ! Fermi energy as reported by CASTEP

  real(kind=dp), allocatable, public, save :: num_electrons(:) ! Holds up-spin and
  ! down-spin

  integer, public, save       :: nbands, nspins
  real(kind=dp), public, save :: electrons_per_state ! 2 for non-spin-P
  ! 1 for spin-P
  logical, public, save       :: spin_polarised

  type, public ::  matrix_weights_array_boundaries
    integer :: norbitals
    integer :: nbands
    integer :: nkpoints
    integer :: nspins
  end type matrix_weights_array_boundaries

  type, public :: orbitals
    integer, allocatable  :: ion_no(:)          ! Unique ion number
    integer, allocatable  :: species_no(:)      ! Unique species number
    integer, allocatable  :: rank_in_species(:) ! Unique ion number within species
    integer, allocatable  :: am_channel(:)      ! The angular momentum Channel (l)
    integer, allocatable  :: shell(:)           ! Principal quantum number (n) !n.b typically only know this for core states
    character(len=10), allocatable :: am_channel_name(:) ! Name of angular momentum channel s,p,d, etc
  end type orbitals

  ! On writing the od2od it became necessary to have some variables moved from
  ! the individual routines into the module.
  integer, allocatable, public, save :: nbands_occ(:, :)
  ! All these are len=80 to be consistent with CASTEP.
  character(len=80), public, save :: omefile_header
  character(len=80), public, save :: domefile_header
  character(len=80), public, save :: pdosfile_header
  character(len=80), public, save :: elnesfile_header
  integer, public, save :: omefile_version
  integer, public, save :: domefile_version
  integer, public, save :: pdosfile_version
  integer, public, save :: elnesfile_version

  type(orbitals), public, save :: pdos_orbital
  real(kind=dp), public, allocatable, save  :: pdos_weights(:, :, :, :)
  real(kind=dp), allocatable :: all_pdos_weights(:, :, :, :)
  type(matrix_weights_array_boundaries), public, save :: pdos_mwab

  type(orbitals), public, save :: elnes_orbital
  type(matrix_weights_array_boundaries), public, save :: elnes_mwab
  real(kind=dp), public, allocatable, save :: all_kpoints(:, :) ! We need this to be available if we're
  ! doing bandgap analysis.
  real(kind=dp), public, allocatable, save :: all_kpoint_weight(:)

  !-------------------------------------------------------------------------!

  private

  !-------------------------------------------------------------------------!
  ! G L O B A L   F U N C T I O N S
  public :: elec_report_parameters
  public :: elec_read_band_energy
  public :: elec_read_band_energy_ordered
  public :: elec_read_band_gradient
  public :: elec_read_optical_mat
  public :: elec_read_elnes_mat
  public :: elec_pdos_read
  public :: elec_pdis_read
  public :: elec_dealloc_elnes
  public :: elec_dealloc_pdos
  public :: elec_dealloc_band_gradient
  public :: elec_dealloc_optical
  public :: elec_elnes_find_channel_names
  public :: elec_elnes_find_channel_numbers

  !-------------------------------------------------------------------------!

contains

  !=========================================================================
  subroutine elec_report_parameters
    !=========================================================================
    ! Report the electronic properties in the calculation
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables nbands,num_electrons,nkpoints,kpoint_grid_dim
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_io, only: stdout
    use od_cell, only: kpoint_grid_dim, nkpoints
    use od_parameters, only: pdis
    implicit none

    write (stdout, *)

    write (stdout, '(1x,a78)') '+----------------------- Electronic Data ------------------------------------+'

    write (stdout, '(1x,a46,i14,a18)') '|  Number of Bands                           :', nbands, "|"
    if (.not. pdis) then
      write (stdout, '(1x,a46,6x,i3,1x,a1,i3,1x,a1,i3,12x,a1)') '|  Grid size                                 :' &
        , kpoint_grid_dim(1), 'x', kpoint_grid_dim(2), 'x', kpoint_grid_dim(3), '|'
    end if
    write (stdout, '(1x,a46,i14,a18)') '|  Number of K-points                        :', nkpoints, "|"

    if (nspins > 1) then
      write (stdout, '(1x,a78)') '|  Spin-Polarised Calculation                :           True                |'
      write (stdout, '(1x,a46,f17.2,a15)') "|  Number of up-spin electrons               :", num_electrons(1), "|"
      write (stdout, '(1x,a46,f17.2,a15)') "|  Number of down-spin electrons             :", num_electrons(2), "|"
    else
      write (stdout, '(1x,a78)') '|  Spin-Polarised Calculation                :           False               |'
      write (stdout, '(1x,a46,f17.2,a15)') "|  Number of electrons                       :", num_electrons(1), "|"

    end if
    write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'

  end subroutine elec_report_parameters

  !=========================================================================
  subroutine elec_read_band_gradient
    !=========================================================================
    ! Read the .cst_ome file in paralell if appropriate. These are the
    ! gradients of the bands at each kpoint.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_gradient,nspins,nbands
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms, only: on_root, my_node_id, num_nodes, root_id,&
         & comms_recv, comms_send, comms_bcast
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_constants, only: bohr2ang, H2eV
    use od_parameters, only: legacy_file_format, iprint, devel_flag
    use od_algorithms, only: algor_dist_array
    implicit none

    integer :: gradient_unit, i, ib, jb, is, ik, inodes, ierr, loop
    character(filename_len) :: gradient_filename
    logical :: exists
    real(kind=dp) :: time0, time1, file_version
    real(kind=dp), parameter :: file_ver = 1.0_dp

    ! Check that we haven't already done this.

    if (allocated(band_gradient)) return

    ! first try to read a velocity file
    if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
      gradient_filename = trim(seedname)//".cst_vel"
    else
      gradient_filename = trim(seedname)//".dome_bin"
    end if
    if (on_root) inquire (file=gradient_filename, exist=exists)
    call comms_bcast(exists, 1)

    if (exists) then  ! good. We are reading from a velocity file

      time0 = io_time()
      if (on_root) then
        if (iprint > 1) write (stdout, '(a)') ' '
        if (iprint > 1) write (stdout, '(a)') ' Reading band gradients from file: '//trim(gradient_filename)
        gradient_unit = io_file_unit()
        if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
          gradient_filename = trim(seedname)//".cst_vel"
          open (unit=gradient_unit, file=gradient_filename, status="old", form='unformatted', err=101)
        else
          gradient_filename = trim(seedname)//".dome_bin"
          open (unit=gradient_unit, file=gradient_filename, status="old", form='unformatted', err=102)
          read (gradient_unit) file_version
          if ((file_version - file_ver) > 0.001_dp) &
            call io_error('Error: Trying to read newer version of dome_bin file. Update optados!')
          read (gradient_unit) domefile_header
          if (iprint > 1) write (stdout, *) trim(domefile_header)
        end if
      end if

      ! Figure out how many kpoint should be on each node
      call algor_dist_array(nkpoints, num_kpoints_on_node)
      allocate (band_gradient(1:nbands, 1:3, 1:num_kpoints_on_node(0), 1:nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: Problem allocating band_gradient in elec_read_band_gradient')

      band_gradient = 0.0_dp
      if (on_root) then
        do inodes = 1, num_nodes - 1
          do ik = 1, num_kpoints_on_node(inodes)
            do is = 1, nspins
              read (gradient_unit) ((band_gradient(ib, i, ik, is), ib=1, nbands), i=1, 3)
            end do
          end do
          call comms_send(band_gradient(1, 1, 1, 1), nbands*3*nspins*num_kpoints_on_node(0), inodes)
        end do
        do ik = 1, num_kpoints_on_node(0)
          do is = 1, nspins
            read (gradient_unit) ((band_gradient(ib, i, ik, is), ib=1, nbands), i=1, 3)
          end do
        end do
      end if

      if (.not. on_root) then
        call comms_recv(band_gradient(1, 1, 1, 1), nbands*3*nspins*num_kpoints_on_node(0), root_id)
      end if

!        write(*,*) "I'm node", my_node_id, "k-pts:", num_kpoints_on_node(my_node_id),"bgarray:", &
!& size(band_gradient), "or:", nbands*3*nspins*num_kpoints_on_node(my_node_id)

      if (on_root) close (unit=gradient_unit)

      ! Convert all band gradients to eV Ang
      band_gradient = band_gradient*bohr2ang*H2eV

      time1 = io_time()
      if (on_root .and. iprint > 1) write (stdout, '(1x,a40,f11.3,a)') 'Time to read band gradients ', time1 - time0, ' (sec)'

    else ! lets try to get the data from the cst_ome file
      allocate (band_gradient(1:nbands, 1:3, 1:num_kpoints_on_node(0), 1:nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: Problem allocating band_gradient (b) in elec_read_band_gradient')

      if (allocated(optical_mat)) then
        do loop = 1, nbands
          band_gradient(loop, :, :, :) = real(optical_mat(loop, loop, :, :, :), dp)
        end do
      else
        call elec_read_optical_mat
        do loop = 1, nbands
          band_gradient(loop, :, :, :) = real(optical_mat(loop, loop, :, :, :), dp)
        end do
        call elec_dealloc_optical  ! given that it is a large matrix we'll let it go
        ! potentially that means reading it twice...
      end if
    end if

    return

101 call io_error('Error: Problem opening cst_vel file in read_band_gradient')
102 call io_error('Error: Problem opening dome_bin file in read_band_gradient')

  end subroutine elec_read_band_gradient

  !=========================================================================
  subroutine elec_read_optical_mat
    !=========================================================================
    ! Read the .cst_ome file in paralell if appropriate. These are the
    ! gradients of the bands at each kpoint.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_gradient,nspins,nbands
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms, only: on_root, my_node_id, num_nodes, root_id,&
         & comms_recv, comms_send
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_constants, only: bohr2ang, H2eV
    use od_parameters, only: legacy_file_format, iprint, devel_flag
    use od_algorithms, only: algor_dist_array
    implicit none

    integer :: gradient_unit, i, ib, jb, is, ik, inodes, ierr
    character(filename_len) :: gradient_filename
    real(kind=dp) :: time0, time1, file_version
    real(kind=dp), parameter :: file_ver = 1.0_dp

    ! Check that we haven't already done this.

    if (allocated(optical_mat)) return

    time0 = io_time()
    if (on_root) then
      gradient_unit = io_file_unit()
      if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
        gradient_filename = trim(seedname)//".cst_ome"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading optical matrix elements from file: '//trim(gradient_filename)
        open (unit=gradient_unit, file=gradient_filename, status="old", form='unformatted', err=101)
      else
        gradient_filename = trim(seedname)//".ome_bin"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading optical matrix elements from file: '//trim(gradient_filename)
        open (unit=gradient_unit, file=gradient_filename, status="old", form='unformatted', err=102)
        read (gradient_unit) file_version
        if ((file_version - file_ver) > 0.001_dp) &
          call io_error('Error: Trying to read newer version of ome_bin file. Update optados!')
        read (gradient_unit) omefile_header
        if (iprint > 1) write (stdout, '(1x,a)') trim(omefile_header)
      end if
    end if

    ! Figure out how many kpoints should be on each node
    call algor_dist_array(nkpoints, num_kpoints_on_node)
    allocate (optical_mat(1:nbands, 1:nbands, 1:3, 1:num_kpoints_on_node(0), 1:nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating optical_mat in elec_read_optical_mat')

    if (legacy_file_format) then

      if (on_root) then
        do inodes = 1, num_nodes - 1
          do ik = 1, num_kpoints_on_node(inodes)
            do is = 1, nspins
              do i = 1, 3
                do jb = 1, nbands
                  do ib = 1, nbands
                    ! Read in units of Ha Bohr^2 / Ang
                    read (gradient_unit) optical_mat(ib, jb, i, ik, is)
                  end do
                end do
              end do
            end do
          end do
          call comms_send(optical_mat(1, 1, 1, 1, 1), nbands*nbands*3*nspins*num_kpoints_on_node(0), inodes)
        end do

        do ik = 1, num_kpoints_on_node(0)
          do is = 1, nspins
            do i = 1, 3
              do jb = 1, nbands
                do ib = 1, nbands
                  ! Read in units of Ha Bohr^2 / Ang
                  read (gradient_unit) optical_mat(ib, jb, i, ik, is)
                end do
              end do
            end do
          end do
        end do
      end if

    else ! sane file format

      if (on_root) then
        do inodes = 1, num_nodes - 1
          do ik = 1, num_kpoints_on_node(inodes)
            do is = 1, nspins
              read (gradient_unit) (((optical_mat(ib, jb, i, ik, is), ib=1, nbands) &
                                     , jb=1, nbands), i=1, 3)
            end do
          end do
          call comms_send(optical_mat(1, 1, 1, 1, 1), nbands*nbands*3*nspins*num_kpoints_on_node(0), inodes)
        end do
        do ik = 1, num_kpoints_on_node(0)
          do is = 1, nspins
            read (gradient_unit) (((optical_mat(ib, jb, i, ik, is), ib=1, nbands), jb=1, nbands), i=1, 3)
          end do
        end do
      end if
    end if

    if (.not. on_root) then
      call comms_recv(optical_mat(1, 1, 1, 1, 1), nbands*nbands*3*nspins*num_kpoints_on_node(0), root_id)
    end if

    if (on_root) close (unit=gradient_unit)

    ! Convert all band gradients to eV Ang
    if (legacy_file_format) then
      optical_mat = optical_mat*bohr2ang*bohr2ang*H2eV
    else
      optical_mat = optical_mat*bohr2ang*H2eV
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to read Optical Matrix Elements                   &
           &      ', time1 - time0, ' (sec) +'
    end if

    return

101 call io_error('Error: Problem opening cst_ome file in read_band_optical_mat')
102 call io_error('Error: Problem opening ome_bin file in read_band_optical_mat')

  end subroutine elec_read_optical_mat

  !=========================================================================
  subroutine elec_read_band_energy !(band_energy,kpoint_r,kpoint_weight)
    !=========================================================================
    ! Read the .bands file in the kpoint list, kpoint weights and band energies
    ! also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_energy, efermi_castep, num_electrons
    ! spin_polarised, electrons_per_state, nspins,nbands
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_constants, only: H2eV
    use od_cell, only: nkpoints, kpoint_r, kpoint_weight, cell_find_MP_grid,&
         & real_lattice, kpoint_grid_dim, num_kpoints_on_node
    use od_comms, only: comms_bcast, comms_send, comms_recv, num_nodes, my_node_id,&
         & on_root, root_id
    use od_io, only: io_file_unit, seedname, filename_len, stdout, io_time,&
         & io_error
    use od_algorithms, only: algor_dist_array
    use od_parameters, only: iprint, compute_band_gap, kpoint_mp_grid

    implicit none

    integer :: inodes, ik, is, ib, band_unit, iall_kpoints, i
    integer :: dum_i1, ierr, str_pos
    character(filename_len) :: band_filename
    character(len=80) :: dummy
    real(kind=dp) :: time0, time1

    time0 = io_time()

    ! Check that we haven't already read in the energies
    if (allocated(band_energy)) return

    !Open the bands file
    band_unit = io_file_unit()
    band_filename = trim(seedname)//".bands"

    ! Read the header from the bands file
    if (on_root) then
      open (unit=band_unit, file=band_filename, status="old", form='formatted', err=100)
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'k-points')
      read (dummy(str_pos + 8:), *) nkpoints
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'components')
      read (dummy(str_pos + 10:), *) nspins
      read (band_unit, '(a)') dummy

      allocate (num_electrons(nspins), stat=ierr)
      if (ierr /= 0) stop " Error : cannot allocate num_electrons"
      str_pos = index(dummy, 'electrons')
      read (dummy(str_pos + 10:), *) num_electrons(:)
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'eigenvalues')
      read (dummy(str_pos + 11:), *) nbands
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'units)')
      read (dummy(str_pos + 6:), '(f12.4)') efermi_castep
      read (band_unit, '(a)') dummy
      read (band_unit, *) real_lattice(:, 1)
      read (band_unit, *) real_lattice(:, 2)
      read (band_unit, *) real_lattice(:, 3)
    end if

    call comms_bcast(nspins, 1)
    if (.not. on_root) then
      allocate (num_electrons(nspins), stat=ierr)
      if (ierr /= 0) stop " Error : cannot allocate num_electrons"
    end if
    call comms_bcast(num_electrons(1), nspins)
    call comms_bcast(nkpoints, 1)
    call comms_bcast(nbands, 1)
    call comms_bcast(efermi_castep, 1)
    !
    call algor_dist_array(nkpoints, num_kpoints_on_node)
    !
    allocate (band_energy(1:nbands, 1:nspins, 1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating band_energy in read_band_energy')
    allocate (kpoint_weight(1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating kpoint_weight in read_band_energy')
    allocate (kpoint_r(1:3, 1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating kpoint_r in read_band_energy')

    if (on_root) then
      allocate (all_kpoints(1:3, nkpoints), stat=ierr)
      if (ierr /= 0) call io_error('Error: Problem allocating all_kpoints in read_band_energy')
      iall_kpoints = 1
      do inodes = 1, num_nodes - 1
        do ik = 1, num_kpoints_on_node(inodes)
          read (band_unit, '(a)') dummy
          str_pos = index(dummy, 'K-point')
          read (dummy(str_pos + 7:), *) dum_i1, kpoint_r(1, ik), kpoint_r(2, ik), kpoint_r(3, ik), kpoint_weight(ik)
          do i = 1, 3
            all_kpoints(i, iall_kpoints) = kpoint_r(i, ik)
          end do
          iall_kpoints = iall_kpoints + 1
          do is = 1, nspins
            read (band_unit, *) dummy
            do ib = 1, nbands
              read (band_unit, *) band_energy(ib, is, ik) !NB spin <-> kpt swapped
            end do
          end do
        end do
        call comms_send(band_energy(1, 1, 1), nbands*nspins*num_kpoints_on_node(0), inodes)
        call comms_send(kpoint_r(1, 1), 3*num_kpoints_on_node(0), inodes)
        call comms_send(kpoint_weight(1), num_kpoints_on_node(0), inodes)
      end do

      do ik = 1, num_kpoints_on_node(0)
        read (band_unit, '(a)') dummy
        str_pos = index(dummy, 'K-point')
        read (dummy(str_pos + 7:), *) dum_i1, kpoint_r(1, ik), kpoint_r(2, ik), kpoint_r(3, ik), kpoint_weight(ik)
        do i = 1, 3
          all_kpoints(i, iall_kpoints) = kpoint_r(i, ik)
        end do
        iall_kpoints = iall_kpoints + 1
        do is = 1, nspins
          read (band_unit, *) dummy
          do ib = 1, nbands
            read (band_unit, *) band_energy(ib, is, ik) !NB spin <-> kpt swapped
          end do
        end do
      end do

      ! Do this here so we can free up the all_kpoints memory, unless we need it to calculate
      ! the kpoints at the band-gap.
      if (kpoint_mp_grid(1) > 0) then
        ! we must have set this manually
        kpoint_grid_dim = kpoint_mp_grid
      else
        call cell_find_MP_grid(all_kpoints, nkpoints, kpoint_grid_dim)
      end if
      if (.not. compute_band_gap) then
        deallocate (all_kpoints, stat=ierr)
        if (ierr /= 0) call io_error('Error: Problem deallocating all_kpoints in read_band_energy')
      end if
    end if

    if (.not. on_root) then
      call comms_recv(band_energy(1, 1, 1), nbands*nspins*num_kpoints_on_node(0), root_id)
      call comms_recv(kpoint_r(1, 1), 3*num_kpoints_on_node(0), root_id)
      call comms_recv(kpoint_weight(1), num_kpoints_on_node(0), root_id)
    end if

    if (on_root) close (unit=band_unit)

    band_energy = band_energy*H2eV
    efermi_castep = efermi_castep*H2eV

    ! Things that follow
    if (nspins .lt. 2) then
      spin_polarised = .false.
      electrons_per_state = 2.0_dp
    else
      spin_polarised = .true.
      electrons_per_state = 1.0_dp
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) write (stdout, '(1x,a40,f11.3,a)') 'Time to read band energies  ', time1 - time0, ' (sec)'

    return
100 call io_error('Error: Problem opening bands file in read_band_energy')
  end subroutine elec_read_band_energy

  !=========================================================================
  subroutine elec_read_band_energy_ordered !(band_energy,kpoint_r,kpoint_weight)
    !=========================================================================
    ! Read the .bands file in the kpoint list, kpoint weights and band energies
    ! also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
    ! in the correct order.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_energy, efermi_castep, num_electrons
    ! spin_polarised, electrons_per_state, nspins,nbands
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_constants, only: H2eV
    use od_cell, only: nkpoints, kpoint_r, kpoint_weight, cell_find_MP_grid,&
         & real_lattice, kpoint_grid_dim, num_kpoints_on_node
    use od_comms, only: comms_bcast, comms_send, comms_recv, num_nodes, my_node_id,&
         & on_root, root_id
    use od_io, only: io_file_unit, seedname, filename_len, stdout, io_time,&
         & io_error
    use od_algorithms, only: algor_dist_array
    use od_parameters, only: iprint, compute_band_gap, kpoint_mp_grid, pdis

    implicit none

    integer :: inodes, ik, is, ib, band_unit, iall_kpoints, i
    integer :: ik_bandfile, ierr, str_pos
    character(filename_len) :: band_filename
    character(len=80) :: dummy
    real(kind=dp) :: time0, time1
    real(kind=dp), allocatable :: all_band_energy(:, :, :)
    real(kind=dp) :: dummy_kpt(4)

    time0 = io_time()

    ! Check that we haven't already read in the energies
    if (allocated(band_energy)) return

    !Open the bands file
    band_unit = io_file_unit()
    band_filename = trim(seedname)//".bands"

    ! Read the header from the bands file
    if (on_root) then
      open (unit=band_unit, file=band_filename, status="old", form='formatted', err=100)
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'k-points')
      read (dummy(str_pos + 8:), *) nkpoints
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'components')
      read (dummy(str_pos + 10:), *) nspins
      read (band_unit, '(a)') dummy

      allocate (num_electrons(nspins), stat=ierr)
      if (ierr /= 0) stop " Error : cannot allocate num_electrons"
      str_pos = index(dummy, 'electrons')
      read (dummy(str_pos + 10:), *) num_electrons(:)
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'eigenvalues')
      read (dummy(str_pos + 11:), *) nbands
      read (band_unit, '(a)') dummy
      str_pos = index(dummy, 'units)')
      read (dummy(str_pos + 6:), '(f12.4)') efermi_castep
      read (band_unit, '(a)') dummy
      read (band_unit, *) real_lattice(:, 1)
      read (band_unit, *) real_lattice(:, 2)
      read (band_unit, *) real_lattice(:, 3)
    end if

    call comms_bcast(nspins, 1)
    if (.not. on_root) then
      allocate (num_electrons(nspins), stat=ierr)
      if (ierr /= 0) stop " Error : cannot allocate num_electrons"
    end if
    call comms_bcast(num_electrons(1), nspins)
    call comms_bcast(nkpoints, 1)
    call comms_bcast(nbands, 1)
    call comms_bcast(efermi_castep, 1)
    !
    call algor_dist_array(nkpoints, num_kpoints_on_node)
    !
    allocate (band_energy(1:nbands, 1:nspins, 1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating band_energy in read_band_energy')
    allocate (all_band_energy(1:nbands, 1:nspins, 1:nkpoints), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating all_band_energy in read_band_energy')
    allocate (kpoint_weight(1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating kpoint_weight in read_band_energy')
    allocate (kpoint_r(1:3, 1:num_kpoints_on_node(0)), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating kpoint_r in read_band_energy')

    if (on_root) then
      allocate (all_kpoints(1:3, nkpoints), stat=ierr)
      allocate (all_kpoint_weight(nkpoints), stat=ierr)
      if (ierr /= 0) call io_error('Error: Problem allocating all_kpoints in read_band_energy')
      do ik = 1, nkpoints
        read (band_unit, '(a)') dummy
        str_pos = index(dummy, 'K-point')
        dummy_kpt = 0
        read (dummy(str_pos + 7:), *) ik_bandfile, dummy_kpt(1), dummy_kpt(2), dummy_kpt(4)
        do i = 1, 3
          all_kpoints(i, ik_bandfile) = dummy_kpt(i)
        end do
        all_kpoint_weight(ik_bandfile) = dummy_kpt(4)
        do is = 1, nspins
          read (band_unit, *) dummy
          do ib = 1, nbands
            read (band_unit, *) all_band_energy(ib, is, ik_bandfile) !NB spin <-> kpt swapped
          end do
        end do
      end do

      ! split bands across node-level arrays
      iall_kpoints = 0
      do inodes = 0, num_nodes - 1
        do ik = 1, num_kpoints_on_node(inodes)
          do i = 1, 3
            kpoint_r(i, ik) = all_kpoints(i, ik + iall_kpoints)
          end do
          kpoint_weight(ik) = all_kpoint_weight(ik + iall_kpoints)
          do is = 1, nspins
            do ib = 1, nbands
              band_energy(ib, is, ik) = all_band_energy(ib, is, ik + iall_kpoints) !NB spin <-> kpt swapped
            end do
          end do
        end do
        iall_kpoints = iall_kpoints + num_kpoints_on_node(inodes)
        ! distribute bands across kpoints
        if (inodes /= 0) then
          call comms_send(band_energy(1, 1, 1), nbands*nspins*num_kpoints_on_node(0), inodes)
          call comms_send(kpoint_r(1, 1), 3*num_kpoints_on_node(0), inodes)
          call comms_send(kpoint_weight(1), num_kpoints_on_node(0), inodes)
        end if
      end do

      ! Do this here so we can free up the all_kpoints memory, unless we need it to calculate
      ! the kpoints at the band-gap or do a pdispersion
      if (kpoint_mp_grid(1) > 0) then
        ! we must have set this manually
        kpoint_grid_dim = kpoint_mp_grid
      else
        call cell_find_MP_grid(all_kpoints, nkpoints, kpoint_grid_dim)
      end if

      if ((.not. compute_band_gap) .and. (.not. pdis)) then
        deallocate (all_kpoints, stat=ierr)
        if (ierr /= 0) call io_error('Error: Problem deallocating all_kpoints in read_band_energy')
        deallocate (all_band_energy, stat=ierr)
        if (ierr /= 0) call io_error('Error: Problem deallocating all_band_energy in read_band_energy')
        deallocate (all_kpoint_weight, stat=ierr)
        if (ierr /= 0) call io_error('Error: Problem deallocating all_kpoint_weight in read_band_energy')
      end if

    end if

    if (.not. on_root) then
      call comms_recv(band_energy(1, 1, 1), nbands*nspins*num_kpoints_on_node(0), root_id)
      call comms_recv(kpoint_r(1, 1), 3*num_kpoints_on_node(0), root_id)
      call comms_recv(kpoint_weight(1), num_kpoints_on_node(0), root_id)
    end if

    if (on_root) close (unit=band_unit)

    band_energy = band_energy*H2eV
    efermi_castep = efermi_castep*H2eV

    ! Things that follow
    if (nspins .lt. 2) then
      spin_polarised = .false.
      electrons_per_state = 2.0_dp
    else
      spin_polarised = .true.
      electrons_per_state = 1.0_dp
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) write (stdout, '(1x,a40,f11.3,a)') 'Time to read band energies  ', time1 - time0, ' (sec)'

    return
100 call io_error('Error: Problem opening bands file in read_band_energy')
  end subroutine elec_read_band_energy_ordered

  !=========================================================================
  subroutine elec_read_elnes_mat
    !=========================================================================
    ! Read the .bands file in the kpoint list, kpoint weights and band energies
    ! also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_energy, efermi_castep, num_electrons
    ! spin_polarised, electrons_per_state, nspins,nbands
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: comms_bcast, comms_send, comms_recv, num_nodes, my_node_id,&
         & on_root, root_id
    use od_io, only: io_file_unit, seedname, filename_len, io_time,&
         & io_error, stdout
    use od_parameters, only: legacy_file_format, devel_flag, iprint

    implicit none

    integer :: inodes, ik, ns, nb, indx
    integer :: ierr, elnes_unit, orb, loop
    character(filename_len) :: elnes_filename
    real(kind=dp) :: time0, time1, file_version
    real(kind=dp), parameter :: file_ver = 1.0_dp

    time0 = io_time()

    ! Check that we haven't already read in the energies
    if (allocated(elnes_mat)) return

    !Open the elnes sfile
    if (on_root) then
      elnes_unit = io_file_unit()
      if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
        elnes_filename = trim(seedname)//".eels_mat"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading elnes matrix elements from file: '//trim(elnes_filename)
        open (unit=elnes_unit, file=elnes_filename, form='unformatted', err=100, status='old')
      else
        elnes_filename = trim(seedname)//".elnes_bin"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading elnes matrix elements from file: '//trim(elnes_filename)
        open (unit=elnes_unit, file=elnes_filename, status="old", form='unformatted', err=102)
        read (elnes_unit) file_version
        if ((file_version - file_ver) > 0.001_dp) &
          call io_error('Error: Trying to read newer version of elnes_bin file. Update optados!')
        read (elnes_unit) elnesfile_header
        if (iprint > 1) write (stdout, '(1x,a)') trim(elnesfile_header)
      end if

      read (elnes_unit) elnes_mwab%norbitals
      read (elnes_unit) elnes_mwab%nbands
      read (elnes_unit) elnes_mwab%nkpoints
      read (elnes_unit) elnes_mwab%nspins

      ! check these agree with band data?

      allocate (elnes_orbital%ion_no(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%ion_no')
      allocate (elnes_orbital%species_no(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%species_no')
      allocate (elnes_orbital%rank_in_species(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbitall%rank_in_species')
      allocate (elnes_orbital%shell(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbitall%shell')
      allocate (elnes_orbital%am_channel(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%am_channel')
      allocate (elnes_orbital%am_channel_name(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%am_channel_name')

      read (elnes_unit) elnes_orbital%species_no(1:elnes_mwab%norbitals)
      read (elnes_unit) elnes_orbital%rank_in_species(1:elnes_mwab%norbitals)
      read (elnes_unit) elnes_orbital%shell(1:elnes_mwab%norbitals)
      read (elnes_unit) elnes_orbital%am_channel(1:elnes_mwab%norbitals)

    end if
    call comms_bcast(elnes_mwab%norbitals, 1)
    call comms_bcast(elnes_mwab%nbands, 1)
    call comms_bcast(elnes_mwab%nkpoints, 1)
    call comms_bcast(elnes_mwab%nspins, 1)
    if (.not. on_root) then
      allocate (elnes_orbital%ion_no(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%ion_no')
      allocate (elnes_orbital%species_no(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%species_no')
      allocate (elnes_orbital%rank_in_species(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbitall%rank_in_species')
      allocate (elnes_orbital%shell(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbitall%shell')
      allocate (elnes_orbital%am_channel(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%am_channel')
      allocate (elnes_orbital%am_channel_name(elnes_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(' Error : elec_read_elnes_mat cannot allocate elnes_orbital%am_channel_name')
    end if
    call comms_bcast(elnes_orbital%species_no(1), elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital%rank_in_species(1), elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital%shell(1), elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital%am_channel(1), elnes_mwab%norbitals)

    ! assume same data distribution as bands

    allocate (elnes_mat(1:elnes_mwab%norbitals, 1:elnes_mwab%nbands, 1:3, &
                        1:num_kpoints_on_node(0), 1:elnes_mwab%nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating elnes_mat in elec_read_elnes_mat')

    if (on_root) then
      if (legacy_file_format) then
        do inodes = 1, num_nodes - 1
          do ik = 1, num_kpoints_on_node(inodes)
            do ns = 1, elnes_mwab%nspins
              do orb = 1, elnes_mwab%norbitals
                do nb = 1, elnes_mwab%nbands
                  read (elnes_unit) (elnes_mat(orb, nb, indx, ik, ns), indx=1, 3)
                end do
              end do
            end do
          end do
          call comms_send(elnes_mat(1, 1, 1, 1, 1), elnes_mwab%norbitals*elnes_mwab%nbands*3* &
                          nspins*num_kpoints_on_node(0), inodes)
        end do

        do ik = 1, num_kpoints_on_node(0)
          do ns = 1, elnes_mwab%nspins
            do orb = 1, elnes_mwab%norbitals
              do nb = 1, elnes_mwab%nbands
                read (elnes_unit) (elnes_mat(orb, nb, indx, ik, ns), indx=1, 3)
              end do
            end do
          end do
        end do

      else ! sane format
        do inodes = 1, num_nodes - 1
          do ik = 1, num_kpoints_on_node(inodes)
            do ns = 1, elnes_mwab%nspins
              read (elnes_unit) (((elnes_mat(orb, nb, indx, ik, ns), orb=1, elnes_mwab%norbitals), &
                                  nb=1, elnes_mwab%nbands), indx=1, 3)
            end do
          end do
          call comms_send(elnes_mat(1, 1, 1, 1, 1), elnes_mwab%norbitals*elnes_mwab%nbands*3*nspins* &
                          num_kpoints_on_node(0), inodes)
        end do

        do ik = 1, num_kpoints_on_node(0)
          do ns = 1, elnes_mwab%nspins
            read (elnes_unit) (((elnes_mat(orb, nb, indx, ik, ns), orb=1, elnes_mwab%norbitals), &
                                nb=1, elnes_mwab%nbands), indx=1, 3)
          end do
        end do
      end if
    end if

    if (.not. on_root) then
      call comms_recv(elnes_mat(1, 1, 1, 1, 1), elnes_mwab%norbitals*elnes_mwab%nbands*3*nspins* &
                      num_kpoints_on_node(0), root_id)
    end if

    if (on_root) close (elnes_unit)

    call elec_elnes_find_channel_names()

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to read Elnes Matrix Elements                     &
           &      ', time1 - time0, ' (sec) +'
    end if

    return

100 call io_error('Error: Problem opening elnes file in elec_read_elnes_mat')
102 call io_error('Error: Problem opening elnes_bin file in elec_read_elnes_mat')

  end subroutine elec_read_elnes_mat

  !=========================================================================
  subroutine elec_elnes_find_channel_names
    !=========================================================================
    !
    ! fill in some extra indexing data
    ! Moved from within elec_read_elnes_mat when I made od2od
    ! AJM 5/12/2019
    use od_io, only: io_error
    implicit none

    integer :: loop

    !  elnes_mwab is a module variable so don't declare.

    do loop = 1, elnes_mwab%norbitals
      if (elnes_orbital%am_channel(loop) == 1) then
        elnes_orbital%am_channel(loop) = 0
        elnes_orbital%am_channel_name(loop) = 's'
      elseif (elnes_orbital%am_channel(loop) == 2) then
        elnes_orbital%am_channel(loop) = 1
        elnes_orbital%am_channel_name(loop) = 'px'
      elseif (elnes_orbital%am_channel(loop) == 3) then
        elnes_orbital%am_channel(loop) = 1
        elnes_orbital%am_channel_name(loop) = 'py'
      elseif (elnes_orbital%am_channel(loop) == 4) then
        elnes_orbital%am_channel(loop) = 1
        elnes_orbital%am_channel_name(loop) = 'pz'
      elseif (elnes_orbital%am_channel(loop) == 5) then
        elnes_orbital%am_channel(loop) = 2
        elnes_orbital%am_channel_name(loop) = 'dzz'
      elseif (elnes_orbital%am_channel(loop) == 6) then
        elnes_orbital%am_channel(loop) = 2
        elnes_orbital%am_channel_name(loop) = 'dzy'
      elseif (elnes_orbital%am_channel(loop) == 7) then
        elnes_orbital%am_channel(loop) = 2
        elnes_orbital%am_channel_name(loop) = 'dzx'
      elseif (elnes_orbital%am_channel(loop) == 8) then
        elnes_orbital%am_channel(loop) = 2
        elnes_orbital%am_channel_name(loop) = 'dxx-yy'
      elseif (elnes_orbital%am_channel(loop) == 9) then
        elnes_orbital%am_channel(loop) = 2
        elnes_orbital%am_channel_name(loop) = 'dxy'
      elseif (elnes_orbital%am_channel(loop) == 10) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fxxx'
      elseif (elnes_orbital%am_channel(loop) == 11) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fyyy'
      elseif (elnes_orbital%am_channel(loop) == 12) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fzzz'
      elseif (elnes_orbital%am_channel(loop) == 13) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fxyz'
      elseif (elnes_orbital%am_channel(loop) == 14) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fz(xx-yy)'
      elseif (elnes_orbital%am_channel(loop) == 15) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fy(zz-xx)'
      elseif (elnes_orbital%am_channel(loop) == 16) then
        elnes_orbital%am_channel(loop) = 3
        elnes_orbital%am_channel_name(loop) = 'fx(yy-zz)'
      else
        call io_error(' Error : unknown angular momentum state in elec_elnes_find_channel_names')
      end if
    end do

  end subroutine elec_elnes_find_channel_names

  !=========================================================================
  subroutine elec_elnes_find_channel_numbers
    !=========================================================================
    !
    ! The elnes_bin has channel numbers 1-16 internally optados thinks about
    ! channel names. So we need to be able to go back and forth.
    !
    ! CASTEP (hence the bin file) and OptaDOS think about am_channel numbers
    ! differently. To keep consistent we convert to CASTEP's numbering scheme
    ! before we write out.
    ! AJM 5/12/2019
    use od_io, only: io_error
    implicit none

    integer :: loop

    do loop = 1, elnes_mwab%norbitals
      selectcase (trim(elnes_orbital%am_channel_name(loop)))
      case ('s')
        elnes_orbital%am_channel(loop) = 1
      case ('px')
        elnes_orbital%am_channel(loop) = 2
      case ('py')
        elnes_orbital%am_channel(loop) = 3
      case ('pz')
        elnes_orbital%am_channel(loop) = 4
      case ('dzz')
        elnes_orbital%am_channel(loop) = 5
      case ('dzy')
        elnes_orbital%am_channel(loop) = 6
      case ('dzx')
        elnes_orbital%am_channel(loop) = 7
      case ('dxx-yy')
        elnes_orbital%am_channel(loop) = 8
      case ('dxy')
        elnes_orbital%am_channel(loop) = 9
      case ('fxxx')
        elnes_orbital%am_channel(loop) = 10
      case ('fyyy')
        elnes_orbital%am_channel(loop) = 11
      case ('fzzz')
        elnes_orbital%am_channel(loop) = 12
      case ('fxyz')
        elnes_orbital%am_channel(loop) = 13
      case ('fz(xx-yy)')
        elnes_orbital%am_channel(loop) = 14
      case ('fy(zz-xx)')
        elnes_orbital%am_channel(loop) = 15
      case ('fx(yy-zz)')
        elnes_orbital%am_channel(loop) = 16
      case default
        call io_error(' Error : unknown angular momentum state in elec_elnes_find_channel_numbers')
      end select
    end do

  end subroutine elec_elnes_find_channel_numbers

  !=========================================================================
  subroutine elec_pdos_read
    !=========================================================================
    ! Read in the full pdos_weights. Write out any variables that we find on the
    ! way. These will be checked for consistency in the dos module. We can't do it
    ! yet as we haven't read the bands file.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: pw, pdos_weights, pdos_orbital
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms, only: comms_bcast, comms_send, comms_recv, num_nodes, my_node_id,&
        & on_root, root_id
    use od_io, only: stdout, io_file_unit, io_error, seedname, filename_len
    use od_cell, only: num_kpoints_on_node
    use od_parameters, only: legacy_file_format, devel_flag, iprint

    implicit none

    ! Band indices used in the read-in of the pdos
    real(kind=dp)                        :: dummyr1, dummyr2, dummyr3
    integer                              :: dummyi, ib, ik, is, iorbitals
    integer                              :: pdos_in_unit, ierr, inodes
    character(filename_len) :: pdos_filename
    real(kind=dp) :: time0, time1, file_version
    real(kind=dp), parameter :: file_ver = 1.0_dp

    if (allocated(pdos_weights)) return

    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R
    if (on_root) then
      pdos_in_unit = io_file_unit()
      if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
        pdos_filename = trim(seedname)//".pdos_weights"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading pdos weights from file: '//trim(pdos_filename)
        open (unit=pdos_in_unit, file=pdos_filename, form='unformatted', err=100, status='old')
      else
        pdos_filename = trim(seedname)//".pdos_bin"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading pdos weights from file: '//trim(pdos_filename)
        open (unit=pdos_in_unit, file=pdos_filename, status="old", form='unformatted', err=102)
        read (pdos_in_unit) file_version
        if ((file_version - file_ver) > 0.001_dp) &
          call io_error('Error: Trying to read newer version of pdos_bin file. Update optados!')
        read (pdos_in_unit) pdosfile_header
        if (iprint > 1) write (stdout, '(1x,a)') trim(pdosfile_header)
      end if

      read (pdos_in_unit) pdos_mwab%nkpoints
      read (pdos_in_unit) pdos_mwab%nspins
      read (pdos_in_unit) pdos_mwab%norbitals
      read (pdos_in_unit) pdos_mwab%nbands
      if (iprint > 3) then
        write (stdout, *) " +==========================================================================+"
        write (stdout, *) " |                          D E B U G   P D O S                             |"
        write (stdout, *) " +==========================================================================+"
        write (stdout, *) " pdos_mwab%nkpoints : ", pdos_mwab%nkpoints
        write (stdout, *) " pdos_mwab%nspins   : ", pdos_mwab%nspins
        write (stdout, *) " pdos_mwab%norbitals: ", pdos_mwab%norbitals
        write (stdout, *) " pdos_mwab%nbands   : ", pdos_mwab%nbands
      endif

      allocate (pdos_orbital%species_no(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%rank_in_species(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%am_channel(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")

      read (pdos_in_unit) pdos_orbital%species_no(1:pdos_mwab%norbitals)
      read (pdos_in_unit) pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
      read (pdos_in_unit) pdos_orbital%am_channel(1:pdos_mwab%norbitals)
      if (iprint > 3) then
        write (stdout, *)
        write (stdout, *) " +--------------------------------------------------------------------------+"
        write (stdout, *) "                                   pdos_orbital%                           "
        write (stdout, *) "       species_no                 rank_in_species                am_channel"
        write (stdout, *) " +--------------------------------------------------------------------------+"
        do iorbitals = 1, pdos_mwab%norbitals
          write (stdout, '(10x,i5,25x,i5,25x,i5)') pdos_orbital%species_no(iorbitals), pdos_orbital%rank_in_species(iorbitals), &
          & pdos_orbital%am_channel(iorbitals)
        end do
        write (stdout, *) " +==========================================================================+ "
      end if
      !-------------------------------------------------------------------------!
    end if
    call comms_bcast(pdos_mwab%norbitals, 1)
    call comms_bcast(pdos_mwab%nbands, 1)
    call comms_bcast(pdos_mwab%nkpoints, 1)
    call comms_bcast(pdos_mwab%nspins, 1)
    if (.not. on_root) then
      allocate (pdos_orbital%species_no(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%rank_in_species(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%am_channel(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
    end if
    call comms_bcast(pdos_orbital%species_no(1), pdos_mwab%norbitals)
    call comms_bcast(pdos_orbital%rank_in_species(1), pdos_mwab%norbitals)
    call comms_bcast(pdos_orbital%am_channel(1), pdos_mwab%norbitals)

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A
    allocate (nbands_occ(1:num_kpoints_on_node(my_node_id), 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate nbands_occ"
    allocate (pdos_weights(1:pdos_mwab%norbitals, 1:pdos_mwab%nbands, &
                           1:num_kpoints_on_node(0), 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate pdos_weights"

    if (on_root) then
      do inodes = 1, num_nodes - 1
        do ik = 1, num_kpoints_on_node(inodes)
          ! The kpoint number, followed by the kpoint-vector
          read (pdos_in_unit) dummyi, dummyr1, dummyr2, dummyr3
          do is = 1, pdos_mwab%nspins
            read (pdos_in_unit) dummyi ! this is the spin number
            read (pdos_in_unit) nbands_occ(ik, is)
            do ib = 1, nbands_occ(ik, is)
              if (iprint > 3) then
                write (stdout, *) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
                write (stdout, *) ib, ik, is
                write (stdout, *) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "
              end if
              read (pdos_in_unit) pdos_weights(1:pdos_mwab%norbitals, ib, ik, is)
            end do
          end do
        end do
        call comms_send(pdos_weights(1, 1, 1, 1), pdos_mwab%norbitals*pdos_mwab%nbands* &
                        nspins*num_kpoints_on_node(0), inodes)
      end do

      do ik = 1, num_kpoints_on_node(0)
        ! The kpoint number, followed by the kpoint-vector
        read (pdos_in_unit) dummyi, dummyr1, dummyr2, dummyr3
        do is = 1, pdos_mwab%nspins
          read (pdos_in_unit) dummyi ! this is the spin number
          read (pdos_in_unit) nbands_occ(ik, is)
          do ib = 1, nbands_occ(ik, is)
            read (pdos_in_unit) pdos_weights(1:pdos_mwab%norbitals, ib, ik, is)
          end do
        end do
      end do
    end if

    if (.not. on_root) then
      call comms_recv(pdos_weights(1, 1, 1, 1), pdos_mwab%norbitals*pdos_mwab%nbands* &
                      nspins*num_kpoints_on_node(0), root_id)
    end if

    if (on_root) close (pdos_in_unit)

    return

100 call io_error('Error: Problem opening pdos_weights file in elec_pdos_read')
102 call io_error('Error: Problem opening pdos_bin file in elec_pdos_read')

  end subroutine elec_pdos_read

  !=========================================================================
  subroutine elec_pdis_read
    !=========================================================================
    ! Read in the full pdos_weights, in correct kpoint path order.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: pw, pdos_weights, pdos_orbital
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Adapted by M L Evans from elec_pdos_read written by A J Morris Aug 2018
    !=========================================================================
    use od_comms, only: comms_bcast, comms_send, comms_recv, num_nodes, my_node_id,&
        & on_root, root_id
    use od_io, only: stdout, io_file_unit, io_error, seedname, filename_len
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_parameters, only: legacy_file_format, devel_flag, iprint

    implicit none

    ! Band indices used in the read-in of the pdos
    integer, allocatable, dimension(:, :) :: all_nbands_occ
    real(kind=dp)                        :: dummyr1, dummyr2, dummyr3
    integer                              :: dummyk, cachek, dummyi, ib, ik, is, stride, iorb
    integer                              :: pdos_in_unit, ierr, inodes, iall_kpoints
    logical                              :: full_debug_pdos_weights = .false.
    character(filename_len) :: pdos_filename
    character(len=80)       :: header
    real(kind=dp) :: time0, time1, file_version
    real(kind=dp), parameter :: file_ver = 1.0_dp

    if (allocated(pdos_weights)) return

    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R
    !-------------------------------------------------------------------------!
    if (on_root) then
      pdos_in_unit = io_file_unit()
      if (index(devel_flag, 'old_filename') > 0 .or. legacy_file_format) then
        pdos_filename = trim(seedname)//".pdos_weights"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading pdos weights from file: '//trim(pdos_filename)
        open (unit=pdos_in_unit, file=pdos_filename, form='unformatted', err=100, status='old')
      else
        pdos_filename = trim(seedname)//".pdos_bin"
        if (iprint > 1) write (stdout, '(1x,a)') 'Reading pdos weights from file: '//trim(pdos_filename)
        open (unit=pdos_in_unit, file=pdos_filename, status="old", form='unformatted', err=102)
        read (pdos_in_unit) file_version
        if ((file_version - file_ver) > 0.001_dp) &
          call io_error('Error: Trying to read newer version of pdos_bin file. Update optados!')
        read (pdos_in_unit) header
        if (iprint > 1) write (stdout, '(1x,a)') trim(header)
      end if

      read (pdos_in_unit) pdos_mwab%nkpoints
      read (pdos_in_unit) pdos_mwab%nspins
      read (pdos_in_unit) pdos_mwab%norbitals
      read (pdos_in_unit) pdos_mwab%nbands
      if (full_debug_pdos_weights) then
        write (stdout, *) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
        write (stdout, *) "pdos_mwab%nkpoints= ", pdos_mwab%nkpoints
        write (stdout, *) "pdos_mwab%nspins= ", pdos_mwab%nspins
        write (stdout, *) "pdos_mwab%norbitals= ", pdos_mwab%norbitals
        write (stdout, *) "pdos_mwab%nbands= ", pdos_mwab%nbands
        write (stdout, *) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "
      end if

      allocate (pdos_orbital%species_no(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%rank_in_species(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%am_channel(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")

      read (pdos_in_unit) pdos_orbital%species_no(1:pdos_mwab%norbitals)
      read (pdos_in_unit) pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
      read (pdos_in_unit) pdos_orbital%am_channel(1:pdos_mwab%norbitals)
      !-------------------------------------------------------------------------!
      call comms_bcast(pdos_mwab%norbitals, 1)
      call comms_bcast(pdos_mwab%nbands, 1)
      call comms_bcast(pdos_mwab%nkpoints, 1)
      call comms_bcast(pdos_mwab%nspins, 1)
    end if
    if (.not. on_root) then
      allocate (pdos_orbital%species_no(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%rank_in_species(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
      allocate (pdos_orbital%am_channel(pdos_mwab%norbitals), stat=ierr)
      if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
    end if
    call comms_bcast(pdos_orbital%species_no(1), pdos_mwab%norbitals)
    call comms_bcast(pdos_orbital%rank_in_species(1), pdos_mwab%norbitals)
    call comms_bcast(pdos_orbital%am_channel(1), pdos_mwab%norbitals)

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A
    allocate (all_nbands_occ(1:nkpoints, 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate all_nbands_occ"
    allocate (pdos_weights(1:pdos_mwab%norbitals, 1:pdos_mwab%nbands, &
                           1:num_kpoints_on_node(0), 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate pdos_weights"
    allocate (all_pdos_weights(1:pdos_mwab%norbitals, 1:pdos_mwab%nbands, &
                               1:nkpoints, 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate all_pdos_weights"

    if (on_root) then
      ! Read in the k-points in the correct path ordering, not the file ordering
      cachek = 0
      stride = 1
      do ik = 1, nkpoints
        ! The kpoint number, followed by the kpoint-vector
        read (pdos_in_unit) dummyk, dummyr1, dummyr2, dummyr3
        if (ik == 2) then
          stride = dummyk - cachek
        end if

        if ((dummyk - cachek) < stride) then
          if (mod(ik, (nkpoints/stride)) == 1) then
            dummyk = dummyk
          else
            dummyk = cachek + stride
          end if
        end if
        cachek = dummyk

        do is = 1, pdos_mwab%nspins
          read (pdos_in_unit) dummyi ! this is the spin number
          read (pdos_in_unit) all_nbands_occ(dummyk, is)
          do ib = 1, all_nbands_occ(dummyk, is)
            read (pdos_in_unit) all_pdos_weights(1:pdos_mwab%norbitals, ib, dummyk, is)
          end do
        end do
      end do

      close (pdos_in_unit)

    end if

    call comms_bcast(all_pdos_weights(1, 1, 1, 1), size(all_pdos_weights))

    iall_kpoints = 0
    do inodes = 0, my_node_id - 1
      iall_kpoints = iall_kpoints + inodes*num_kpoints_on_node(inodes)
    end do

    do ik = 1, num_kpoints_on_node(my_node_id)
      do is = 1, pdos_mwab%nspins
        do ib = 1, all_nbands_occ(ik + iall_kpoints, is)
          do iorb = 1, pdos_mwab%norbitals
            pdos_weights(iorb, ib, ik, is) = all_pdos_weights(iorb, ib, ik + iall_kpoints, is)
          end do
        end do
      end do
    end do
    deallocate (all_pdos_weights)

    return

100 call io_error('Error: Problem opening pdos_weights file in elec_pdis_read')
102 call io_error('Error: Problem opening pdos_bin file in elec_pdis_read')

  end subroutine elec_pdis_read

  subroutine elec_dealloc_pdos
    use od_io, only: io_error
    implicit none
    integer :: ierr

    if (allocated(pdos_weights)) then
      deallocate (pdos_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating pdos_weights in elec_dealloc_pdos')
    end if

  end subroutine elec_dealloc_pdos

  subroutine elec_dealloc_elnes
    use od_io, only: io_error
    implicit none
    integer :: ierr

    if (allocated(elnes_mat)) then
      deallocate (elnes_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating elnes_mat in elec_dealloc_elnes')
    end if

!    if(allocated(elnes_orbital)) then
!       deallocate(elnes_orbital,stat=ierr)
!       if (ierr/=0) call io_error('Error in deallocating elnes_orbital in elec_dealloc_elnes')
!    end if

  end subroutine elec_dealloc_elnes

  subroutine elec_dealloc_band_gradient
    use od_io, only: io_error
    implicit none
    integer :: ierr

    if (allocated(band_gradient)) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating band_gradient in elec_dealloc_band_gradient')
    end if

  end subroutine elec_dealloc_band_gradient

  subroutine elec_dealloc_optical
    use od_io, only: io_error
    implicit none
    integer :: ierr

    if (allocated(optical_mat)) then
      deallocate (optical_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error in deallocating optical_mat in elec_dealloc_optical')
    end if

  end subroutine elec_dealloc_optical

end module od_electronic
