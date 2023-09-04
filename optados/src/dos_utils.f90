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
! MODULE od_dos  OptaDOS - Density of States
! This is the module that contains all if the DOS routines. It is used through
! the global calculate_dos subroutine
!-------------------------------------------------------------------------------
! Three other global routines are available, dos_merge, doslin and
! doslin_sub_cell_corners, these are currently used by the jdos module and
! routines.
!-------------------------------------------------------------------------------
! Written by: A J Morris Nov - Dec 2010 Modified from LinDOS (CJP+AJM)
!===============================================================================
module od_dos_utils
  use od_constants, only: dp
  use od_electronic, only: matrix_weights_array_boundaries
  implicit none

  !-------------------------------------------------------------------------------
  ! D E R I V E D   P R O T O T Y P E S

  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! P U B L I C   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: dos_adaptive(:, :)
  real(kind=dp), allocatable, public, save :: dos_fixed(:, :)
  real(kind=dp), allocatable, public, save :: dos_linear(:, :)

  real(kind=dp), allocatable, public, save :: intdos_adaptive(:, :)
  real(kind=dp), allocatable, public, save :: intdos_fixed(:, :)
  real(kind=dp), allocatable, public, save :: intdos_linear(:, :)

  real(kind=dp), allocatable, public, save :: E(:)
  real(kind=dp), public, save :: vbm_energy = 0.0_dp
  real(kind=dp), public, save :: cbm_energy = 0.0_dp

  real(kind=dp), public, save :: efermi_fixed
  real(kind=dp), public, save :: efermi_adaptive
  real(kind=dp), public, save :: efermi_linear

  !real(kind=dp), public, save :: efermi_quad
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S
  public :: dos_utils_calculate
  public :: dos_utils_deallocate
  public :: dos_utils_calculate_at_e
  public :: dos_utils_merge          ! Used by od_jdos
  public :: doslin_sub_cell_corners  ! Used by od_jdos
  public :: doslin                   ! Used by od_jdos
  public :: dos_utils_compute_dos_at_efermi
  public :: dos_utils_compute_bandgap
  public :: dos_utils_compute_band_energies
  public :: dos_utils_set_efermi
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! P R I V A T E   V A R I A B L E S
  private ! unless otherwise indicated

  type(matrix_weights_array_boundaries) :: mw

  real(kind=dp), save                   :: delta_bins ! Width of bins
  logical :: calc_weighted_dos
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine dos_utils_calculate(matrix_weights, weighted_dos)
    !===============================================================================
    ! Main routine in dos module, drives the calculation of density of states for
    ! both task : dos and also if it is required elsewhere.
    !-------------------------------------------------------------------------------
    ! Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
    !            weighted_dos   (out)(opt) : Output DOS weigthed by matrix_weights
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
    ! intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
    ! efermi_linear, delta_bins, calc_weighted_dos
    !-------------------------------------------------------------------------------
    ! Modules Used: see below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: One of linear, adaptive or fixed must be .true.
    !-------------------------------------------------------------------------------
    ! Known Worries: (1) If more than one of linear, adaptive or fixed are set it
    ! uses the most complicated method.
    ! (2) It should be possible to pass optioinal arguments to sub programs as
    ! optional argumnets without checking whether they are there or not. g95 will
    ! allow this behaviour. gfotran will not.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !===============================================================================
    use od_io, only: stdout, io_time, io_error
    use od_comms, only: on_root, my_node_id, comms_bcast
    use od_electronic, only: band_gradient, band_energy, efermi, efermi_castep, nspins, &
         & elec_read_band_gradient, unshifted_efermi
    use od_parameters, only: linear, adaptive, fixed, quad, &
         & dos_per_volume, iprint, set_efermi_zero, efermi_choice, iprint, photo, photo_slab_volume
    use od_cell, only: cell_volume, num_kpoints_on_node

    implicit none

    !-------------------------------------------------------------------------------
    ! I N T E R N A L   V A R I A B L E S
    real(kind=dp) :: time0, time1
    real(kind=dp), intent(in), allocatable, optional  :: matrix_weights(:, :, :, :)
    real(kind=dp), intent(out), allocatable, optional  :: weighted_dos(:, :, :) ! bins.spins, orbitals

    !-------------------------------------------------------------------------------

    if (.not. (linear .or. adaptive .or. fixed .or. quad)) call io_error(" DOS: No Broadening Set")

    calc_weighted_dos = .false.
    if (present(matrix_weights)) calc_weighted_dos = .true.

    if (calc_weighted_dos .eqv. .false.) then ! We are called just to provide dos.
      if (allocated(E)) then
        if (on_root .and. iprint > 1) write (stdout, '(1x,a78)') "| Already calculated dos, so returning... &
             &                                   |"
        return  ! The dos has already been calculated previously so just return.
      end if
    end if

    if (calc_weighted_dos) then
      mw%norbitals = size(matrix_weights, 1)
      mw%nbands = size(matrix_weights, 2)
      mw%nkpoints = size(matrix_weights, 3)
      mw%nspins = size(matrix_weights, 4)
    end if

    if (calc_weighted_dos) then
      !       print*,'mw%nkpoints.ne.num_nkpoints_on_node(my_node_id))',mw%nkpoints,nunum_nkpoints_on_node(my_node_id)
      if (mw%nspins .ne. nspins) call io_error("ERROR : DOS :  mw%nspins not equal to nspins.")
      if (mw%nkpoints .ne. num_kpoints_on_node(my_node_id)) &
        call io_error("ERROR : DOS : mw%nkpoints not equal to nkpoints.")
    end if

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if (quad .or. linear .or. adaptive) then
      if (.not. allocated(band_gradient)) call elec_read_band_gradient
    end if
    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   D O S
    ! Now everything is set up, we can perform the dos accumulation in parallel
    time0 = io_time()

    call setup_energy_scale

    if (on_root .and. (iprint > 1)) write (stdout, *)

    if (fixed) then
      if (calc_weighted_dos .and. (.not. adaptive) .and. (.not. linear)) then
        call calculate_dos("f", dos_fixed, intdos_fixed, matrix_weights=matrix_weights, weighted_dos=weighted_dos)
        call dos_utils_merge(dos_fixed, weighted_dos=weighted_dos)
      else
        call calculate_dos("f", dos_fixed, intdos_fixed)
        call dos_utils_merge(dos_fixed)
      end if
      call dos_utils_merge(intdos_fixed)
    end if
    if (adaptive) then
      if (calc_weighted_dos .and. (.not. linear)) then
        call calculate_dos("a", dos_adaptive, intdos_adaptive, matrix_weights=matrix_weights, weighted_dos=weighted_dos)
        call dos_utils_merge(dos_adaptive, weighted_dos=weighted_dos)
      else
        call calculate_dos("a", dos_adaptive, intdos_adaptive)
        call dos_utils_merge(dos_adaptive)
      end if
      call dos_utils_merge(intdos_adaptive)
    end if
    if (linear) then
      if (calc_weighted_dos) then
        call calculate_dos("l", dos_linear, intdos_linear, matrix_weights=matrix_weights, weighted_dos=weighted_dos)
        call dos_utils_merge(dos_linear, weighted_dos=weighted_dos)
      else
        call calculate_dos("l", dos_linear, intdos_linear)
        call dos_utils_merge(dos_linear)
      end if
      call dos_utils_merge(intdos_linear)
    end if

    if (quad) then
      call io_error("quadratic broadening not implemented")
      !if(quad)    call merge_dos(dos_quad)
      !if(quad)    call merge_dos(intdos_quad)
    end if

    !    if(.not.on_root) then
    !       if(allocated(E)) deallocate(E, stat=ierr)
    !       if (ierr/=0) call io_error ("cannot deallocate  E")
    !    endif

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a59,f11.3,a8)') '+ Time to calculate DOS                                     ', &
      &time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)')
    end if
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! F E R M I   E N E R G Y   A N A L Y S I S
    if (efermi_choice == "optados") then
      if (on_root) then
        time0 = io_time()
        write (stdout, '(1x,a78)') '+----------------------------- Fermi Energy Analysis ------------------------+'
        !    write(stdout,'(1x,a1,a46,f8.4,a3,12x,a8)') "|"," Fermi energy from CASTEP : ",efermi_castep," eV","<- EfC |"
        !    write(stdout,'(1x,a71)')  '+---------------------------------------------------------------------+'

        if (fixed) then
          write (stdout, '(1x,a23,54x,a1)') "| From Fixed broadening", "|"
          efermi_fixed = calc_efermi_from_intdos(intdos_fixed)
          write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|", " Fermi energy (Fixed broadening) : ", &
          &efermi_fixed, "eV", "<- EfF |"

          write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        end if
        if (adaptive) then
          write (stdout, '(1x,a26,51x,a1)') "| From Adaptive broadening", "|"
          efermi_adaptive = calc_efermi_from_intdos(intdos_adaptive)
          write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|", " Fermi energy (Adaptive broadening) : " &
            , efermi_adaptive, "eV", "<- EfA |"

          write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        end if
        if (linear) then
          write (stdout, '(1x,a24,53x,a1)') "| From Linear broadening", "|"
          efermi_linear = calc_efermi_from_intdos(intdos_linear)
          write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|", " Fermi energy (Linear broadening) : ", &
            efermi_linear, " eV", "<- EfL |"

          write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        end if
      end if
    end if

    call comms_bcast(efermi_fixed, 1)
    call comms_bcast(efermi_linear, 1)
    call comms_bcast(efermi_adaptive, 1)

    if (on_root) then
      if (dos_per_volume) then
        if (fixed) then
          dos_fixed = dos_fixed/cell_volume
          intdos_fixed = intdos_fixed/cell_volume
        end if
        if (adaptive) then
          dos_adaptive = dos_adaptive/cell_volume
          intdos_adaptive = intdos_adaptive/cell_volume
        end if
        if (linear) then
          dos_linear = dos_linear/cell_volume
          intdos_linear = intdos_linear/cell_volume
        end if
        if (calc_weighted_dos) then
          if (photo) then
            weighted_dos = weighted_dos/photo_slab_volume
          else
            weighted_dos = weighted_dos/cell_volume
          end if
        end if
        ! if(quad) then
        !    dos_quad=dos_quad/cell_volume
        !    intdos_quad=intdos_quad/cell_volume
        ! endif
      end if
    end if

  end subroutine dos_utils_calculate

  !===============================================================================
  subroutine dos_utils_set_efermi
    !===============================================================================
    use od_parameters, only: efermi_choice, efermi_user, fixed,&
         & linear, adaptive, iprint
    use od_electronic, only: efermi_castep, num_electrons, nspins, efermi, &
         & electrons_per_state, band_energy, nbands, efermi_set
    use od_io, only: io_error, stdout
    use od_comms, only: on_root, my_node_id, comms_reduce, comms_bcast
    use od_cell, only: num_kpoints_on_node
    implicit none

    integer :: is, ik, top_occ_band
    real(kind=dp) :: vbm, cbm

    if (on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Setting Fermi Energy  --------------------------+'
    end if

!    if(on_root) write(stdout,'(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
!         &" Fermi energy from file : ",efermi_castep," eV","| <- EfC"

    select case (efermi_choice)
    case ("file")
      if (on_root) write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
           &" Set fermi energy from file : ", efermi_castep, " eV", "  <- EfC"
      efermi = efermi_castep
    case ("user")
      if (on_root) write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
           &" Fermi energy from user : ", efermi_user, " eV", "  <- EfU"
      efermi = efermi_user
    case ("insulator")
      ! Same fermi level for up and down spins. Different number
      ! of electrons for up and down spins.
      ! For an insulator. Hence same number of electrons at each kpoint.
      ! band_energy(ib,is,ik)
      vbm = -huge(vbm)
      cbm = huge(cbm)
      if (on_root .and. iprint > 3) write (stdout, *) vbm, " =vbm : cbm= ", cbm

      ! Go between global VBM and CBM band_energy(ib,is,ik)
      do ik = 1, num_kpoints_on_node(my_node_id)
        do is = 1, nspins
          ! Which is the band below the fermi energy at this spin and kpoint
          top_occ_band = ceiling(num_electrons(is)/electrons_per_state)

          if (band_energy(top_occ_band, is, ik) > vbm) &
               &vbm = band_energy(top_occ_band, is, ik)
          ! If the band_energy array is big enough then there will be occupied states.
          if (num_electrons(is) + 1 .le. nbands) then
            if (band_energy(top_occ_band + 1, is, ik) < cbm) &
                 &cbm = band_energy(top_occ_band + 1, is, ik)
          end if
        end do
        if (on_root .and. iprint > 3) write (stdout, *) vbm, " =vbm : cbm= ", cbm
      end do

      ! Find the globals
      call comms_reduce(cbm, 1, 'MIN')
      call comms_bcast(cbm, 1)
      call comms_reduce(vbm, 1, 'MAX')
      call comms_bcast(vbm, 1)

      ! If we have a CBM then set the efermi halfway between VBM and CBM
      ! If we don't have a CBM then set it 0.5 eV above VBM and hope for
      ! the best.
      if (cbm == huge(cbm)) then
        efermi = vbm + 0.5_dp
      else
        efermi = vbm + 0.5_dp*(cbm - vbm)
      end if

      if (on_root) write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
           &" Fermi energy assuming insulator : ", efermi, " eV", "  <- EfI"

    case ("optados")
      ! So in the case of compare_jdos we pick efermi_adaptive.
      call dos_utils_calculate ! This will return if we already have.
      if (fixed) efermi = efermi_fixed
      if (linear) efermi = efermi_linear
      if (adaptive) efermi = efermi_adaptive
      if (on_root) write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
           &" Fermi energy from DOS : ", efermi, " eV", "<- EfD |"
    case default
      call io_error('Error in dos_utils_set_efermi: unknown efermi choice this is a bug')
    end select
    if (on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a78)')
    end if

    efermi_set = .true.

  end subroutine dos_utils_set_efermi

  !===============================================================================
  subroutine dos_utils_compute_dos_at_efermi
    !===============================================================================
    use od_io, only: stdout, io_time
    use od_comms, only: on_root, comms_bcast
    use od_electronic, only: efermi, nspins
    use od_parameters, only: fixed, linear, adaptive, iprint, compute_band_gap
    implicit none

    real(dp) :: time0, time1
    real(dp) :: dos_at_efermi(1:3, 1:nspins) ! Fix,Adapt,Linear
    integer :: is

    time0 = io_time()

    call dos_utils_calculate_at_e(efermi, dos_at_e=dos_at_efermi)

    if ((iprint > 1) .and. on_root) then
      write (stdout, *)
    end if

    if (on_root) then
      write (stdout, '(1x,a78)') '+----------------------- DOS at Fermi Energy Analysis -----------------------+'
      write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|", " Fermi energy used : ", efermi, "eV", "       |"

      if (fixed) then
        write (stdout, '(1x,a78)') "| From Fixed broadening                                                      |"
        do is = 1, nspins
          write (stdout, '(1x,a1,a20,i1,a25,f8.4,a9,6x,a8)') "|", "Spin Component : ", is,&
               &"  DOS at Fermi Energy : ", dos_at_efermi(1, is), " eln/cell", "<- DEF |"
        end do
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if

      if (adaptive) then
        write (stdout, '(1x,a78)') "| From Adaptive broadening                                                   |"
        do is = 1, nspins
          write (stdout, '(1x,a1,a20,i1,a25,f8.4,a9,6x,a8)') "|", "Spin Component : ", is,&
               &"  DOS at Fermi Energy : ", dos_at_efermi(2, is), " eln/cell", "<- DEA |"
        end do
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if

      if (linear) then
        write (stdout, '(1x,a78)') "| From Linear broadening                                                     |"
        do is = 1, nspins
          write (stdout, '(1x,a1,a20,i1,a25,f8.4,a9,6x,a8)') "|", "Spin Component : ", is,&
               &"  DOS at Fermi Energy : ", dos_at_efermi(3, is), " eln/cell", "<- DEL |"
        end do
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      end if

      time1 = io_time()
      if (iprint > 1) then
        write (stdout, '(1x,a59,f11.3,a8)') &
             '+ Time to calculate DOS at Fermi Energy                  &
             &      ', time1 - time0, ' (sec) +'

      end if
      write (stdout, '(1x,a78)')

      !-------------------------------------------------------------------------------
    end if
    call comms_bcast(compute_band_gap, 1)
  end subroutine dos_utils_compute_dos_at_efermi

  !===============================================================================
  subroutine dos_utils_compute_bandgap
    !===============================================================================
    ! Modified from LINDOS -- AJM 3rd June 2011
    ! Rewritten 31/1/12 AJM
    use od_electronic, only: nspins, nbands, efermi, band_energy, num_electrons, &
         &  all_kpoints
    use od_cell, only: nkpoints, num_kpoints_on_node
    use od_io, only: stdout, io_time, io_error
    use od_parameters, only: iprint
    use od_comms, only: comms_send, comms_recv, on_root, my_node_id, &
         &num_nodes, root_id
    implicit none

    ! Local array containing vbm and cbm at each kpoint and spin.
    real(dp), allocatable :: bandgap(:, :, :)        !nbands=1,2,nspins,nkpoints

    ! Same as bandgap only merged over all nodes to the root-node
    real(dp), allocatable :: global_bandgap(:, :, :)

    ! Avergae bandgap for up and down spins (summed over all kpoints and divided
    ! by num_kpoints)
    real(dp), allocatable :: average_bandgap(:)

    ! Optical bandgap. Smallest vertical distance between two kpoints
    real(dp), allocatable ::  optical_bandgap(:)

    ! The number of kpoints that have the vbm and cbm
    integer, allocatable :: thermal_vbm_multiplicity(:), thermal_cbm_multiplicity(:)

    ! The number of kpoints that have the optical gap
    integer, allocatable :: optical_multiplicity(:)

    ! Timing variables
    real(dp) :: time0, time1

    ! Loop counters
    integer :: ik, is, ib, ierr, inode, kpoints_before_this_node

    ! Temporary variables
    real(dp), allocatable :: kpoint_bandgap(:)
    real(dp) :: thermal_cbm, thermal_vbm, thermal_bandgap, weighted_average
    integer :: thermal_vbm_k, thermal_cbm_k
    logical :: thermal_multiplicity

    time0 = io_time()

    ! Preamble
    if (on_root) then
      if (iprint > 2) then
        write (stdout, *)
        write (stdout, '(1x,a78)') "| Finding an estimate of the maximum bandgap...                              |"
      end if
      write (stdout, '(1x,a78)') '+----------------------------- Bandgap Analysis -----------------------------+'
    end if

    if (.not. allocated(all_kpoints)) call io_error('Error all_kpoints not allocated in dos_utils: compute_bandgap')

    allocate (bandgap(1:2, 1:nspins, 1:num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating bandgap in dos_utils: compute_bandgap')
    bandgap = -200.00_dp

    ! Write the VBM and CBM into a slice of an band_energy array
    ! It is easier to merge these slices onto the root node rather than try
    ! the horrendous book-keeping of which kpoint on which node has which
    ! gap
    do is = 1, nspins
      do ik = 1, num_kpoints_on_node(my_node_id)
        bands_loop: do ib = 1, nbands
          ! If this eignevalue is greater than efermi, then this is the CBM
          ! and the one below it is the VBM
          if (band_energy(ib, is, ik) .gt. efermi) then
            bandgap(1, is, ik) = band_energy(ib - 1, is, ik)
            bandgap(2, is, ik) = band_energy(ib, is, ik)
            exit bands_loop
          end if
        end do bands_loop
      end do !ik
    end do !is

    ! Create a global kpoint array,and set it to something silly in case we need to debug whether
    ! data was actually written to it.
    if (on_root) then
      allocate (global_bandgap(1:2, 1:nspins, 1:nkpoints), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating global_bandgap in dos_utils: compute_bandgap')
      global_bandgap = -100.00_dp
    end if

    ! Pass all the slices around Efermi to the head node, making sure whe get the global
    ! kpoint number. *And* crucially the same kpoint numbers as in the bands file.
    ! Otherwise the kpoint numbers of the VBM and CBM change as different numbers of nodes
    ! are used.
    kpoints_before_this_node = 0
    do inode = 1, (num_nodes - 1)
      if (my_node_id == inode) call comms_send(bandgap(1, 1, 1), 2*is*num_kpoints_on_node(inode), root_id)
      if (on_root) call comms_recv(global_bandgap(1, 1, &
           & kpoints_before_this_node + 1), 2*is*num_kpoints_on_node(inode), inode)
      kpoints_before_this_node = kpoints_before_this_node + num_kpoints_on_node(inode)
    end do
    ! Copy the root node's slice to the global array.
    if (on_root) &
         & global_bandgap(1:2, 1:nspins, kpoints_before_this_node &
         & + 1:kpoints_before_this_node + num_kpoints_on_node(root_id)) &
             & = bandgap(1:2, 1:nspins, 1:num_kpoints_on_node(root_id))

    ! We now have an array with the VBM and CBM of all kpoints, so we can deallocate the local ones.
    deallocate (bandgap, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating bandgap in dos_utils: compute_bandgap')

    ! We've freed up a bit of memory, so now we can allocate all of the output data arrays on the
    ! headnode
    if (on_root) then
      allocate (average_bandgap(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating average_bandgap in dos_utils: compute_bandgap')
      allocate (optical_bandgap(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating optical_bandgap in dos_utils: compute_bandgap')
      allocate (kpoint_bandgap(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating kpoint_bandgap in dos_utils: compute_bandgap')
      allocate (thermal_vbm_multiplicity(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating thermal_vbm_multiplicity in dos_utils: compute_bandgap')
      allocate (thermal_cbm_multiplicity(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocatingthermal_cbm_multiplicity in dos_utils: compute_bandgap')
      allocate (optical_multiplicity(nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating optical_multiplicity in dos_utils: compute_bandgap')
    end if

    ! Now lets look for the various bandgaps by cycling through the global_bandgap array
    if (on_root) then
      ! COMMENTED OUT: This writes out the complete global_bandgap array to fort.7 if you'd like
      ! to take a look
      !  do ik=1,nkpoints
      !     do is=1,nspins
      !        write(7,*) ik,is, global_bandgap(1,is,ik),global_bandgap(2,is,ik)
      !     enddo
      ! enddo
      average_bandgap = 0.0_dp
      thermal_cbm = huge(1.0_dp)
      thermal_vbm = -huge(1.0_dp)
      thermal_vbm_multiplicity = 1
      thermal_cbm_multiplicity = 1
      thermal_vbm_k = -1
      thermal_cbm_k = -1
      optical_multiplicity = 1
      optical_bandgap = huge(1.0_dp)

      do is = 1, nspins
        do ik = 1, nkpoints
          ! Look for the thermal vbm. This is done for both spin components although
          ! ws care about the individual spin multiplicity
          if (global_bandgap(1, is, ik) .ge. thermal_vbm) then
            if (abs(global_bandgap(1, is, ik) - thermal_vbm) < epsilon(thermal_vbm)) then
              ! If this is the same value of vbm as we had before, then
              ! we have multiple maxima, and we need to know that we might not
              ! get the direct/indirect gap right
              thermal_vbm_multiplicity(is) = thermal_vbm_multiplicity(is) + 1
            else
              ! If we haven't had this high a vbm value before, then we take it
              ! and set the multiplicty back to zero
              thermal_vbm = global_bandgap(1, is, ik)
              thermal_vbm_multiplicity(is) = 1
              thermal_vbm_k = ik
            end if
          end if
          ! We search for the CBM in the same way as the VBM above.
          if (global_bandgap(2, is, ik) .le. thermal_cbm) then
            if (abs(global_bandgap(2, is, ik) - thermal_cbm) < epsilon(thermal_cbm)) then
              thermal_cbm_multiplicity(is) = thermal_cbm_multiplicity(is) + 1
            else
              thermal_cbm = global_bandgap(2, is, ik)
              thermal_cbm_multiplicity(is) = 1
              thermal_cbm_k = ik
            end if
          end if

          ! Work out the bandgap for this particular kpoint and spin
          kpoint_bandgap(is) = global_bandgap(2, is, ik) - global_bandgap(1, is, ik)

          ! If this is smaller than the optical gap. Then this is our next iteration
          ! of the optical gap. Worry about mutiplicties in the same way. Although this
          ! is just for into, since we can't have direct/indirect optical gaps.
          if (kpoint_bandgap(is) .le. optical_bandgap(is)) then
            if (abs(kpoint_bandgap(is) - optical_bandgap(is)) < epsilon(optical_bandgap(is))) then
              optical_multiplicity(is) = optical_multiplicity(is) + 1
            else
              optical_bandgap(is) = kpoint_bandgap(is)
              optical_multiplicity(is) = 1
            end if
          end if
          average_bandgap(is) = average_bandgap(is) + (global_bandgap(2, is, ik) - global_bandgap(1, is, ik))
        end do ! nkpoints
      end do ! nspins

      ! We now have enough info to calculate the average bandgap and the thermal bandgap
      average_bandgap = average_bandgap/nkpoints

      vbm_energy = thermal_vbm
      cbm_energy = thermal_cbm
      thermal_bandgap = thermal_cbm - thermal_vbm

      ! Don't want this array hanging around a moment longer than we need it.
      deallocate (global_bandgap, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating global_bandgap in dos_utils: compute_bandgap')

      ! If there is more than one VBM and CBM let's flag it up.
      ! We wouldn't like to comment on the nature of the gap at this point
      ! Unless we did some more work. (Which should have been done before we deallocated
      ! global_bandgap!)
      thermal_multiplicity = .false.
      do is = 1, nspins
        if ((thermal_vbm_multiplicity(is) .ne. 1) .or. (thermal_cbm_multiplicity(is) .ne. 1)) thermal_multiplicity = .true.
      end do

      ! Report the thermal gap multiplicity
      write (stdout, '(1x,1a,a50, 19x, a8)') "|", "Number of kpoints at       VBM       CBM", "       |"
      do is = 1, nspins
        write (stdout, '(1x,a1,a25,1x,i3,1x,a3,1x,i5,5x,i5,20x,a8)') "|", " Spin :", is, " : ", &
             &thermal_vbm_multiplicity(is), thermal_cbm_multiplicity(is), "       |"
      end do

      ! Write out the thermal gap info
      write (stdout, '(1x,a1,a32,f15.10,1x,a3,18x,a8)') "|", "Thermal Bandgap :", thermal_bandgap, "eV", "<- TBg |"

      if (.not. thermal_multiplicity) then
        if (thermal_vbm_k == thermal_cbm_k) then
          write (stdout, '(1x,a1,a32,1x,f10.5,1x,f10.5,1x,f10.5,4x,a8)') "|", "At kpoint :", &
            all_kpoints(1:3, thermal_vbm_k), "       |"
          write (stdout, '(1x,a78)') '|             ==> Direct Gap                                                  |'
        else
          write (stdout, '(1x,a1,a32,1x,f10.5,1x,f10.5,1x,f10.5,4x,a8)') "|", "Between VBM kpoint :", &
            all_kpoints(1:3, thermal_vbm_k), "       |"
          write (stdout, '(1x,a1,a32,1x,f10.5,1x,f10.5,1x,f10.5,4x,a8)') "|", "and CBM kpoint:", &
            all_kpoints(1:3, thermal_cbm_k), "       |"
          write (stdout, '(1x,a78)') '|             ==> Indirect Gap                                               |'
        end if
      else ! thermal_mutiplicty=.true.
        write (stdout, '(1x,a78)') '|          ==> Multiple Band Minima and Maxima -- Gap unknown                |'
      end if

      ! We allocated this in elec_read_band_energy but kept it if compute_band_gap=T
      deallocate (all_kpoints, stat=ierr)
      if (ierr /= 0) call io_error('Error: Problem deallocating all_kpoints in read_band_energy')

      ! Write out the Optical gap indo
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a1,a45,a32)') "|", "Optical Bandgap  ", "|"
      do is = 1, nspins
        write (stdout, '(1x,a1,a25,1x,i3,1x,a3,1x,f15.10,1x,a3,16x,a8)') "|", " Spin :", is, " : ", optical_bandgap(is),&
             &"eV", "<- OBg |"
      end do
      write (stdout, '(1x,1a,a50, 19x, a8)') "|", "Number of kpoints with this gap         ", "       |"
      ! The multiplicity info here is just for reference.
      do is = 1, nspins
        write (stdout, '(1x,a1,a25,1x,i3,1x,a3,6x,i5,25x,a8)') "|", " Spin :", is, " : ", optical_multiplicity(is), "       |"
      end do

      ! Write out the average bandgap
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a1,a45,a32)') "|", "Average Bandgap  ", "|"
      do is = 1, nspins
        write (stdout, '(1x,a1,a25,1x,i3,1x,a3,1x,f15.10,1x,a3,16x,a8)') "|", " Spin :", is, " : ", &
             &average_bandgap(is), "eV", "<- ABg |"
      end do
      ! If we have more then one spin, then we need some way to combine the up and down spin bandgaps
      ! At Richard Needs' suggestion we use the weighted sum.
      if (nspins > 1) then
        weighted_average = (average_bandgap(1)*num_electrons(1) + average_bandgap(2)*num_electrons(2)) &
             & /(num_electrons(1) + num_electrons(2))
        write (stdout, '(1x,a1,a33,1x,f15.10,1x,a3,16x,a8)') "|", " Weighted Average : ", weighted_average, "eV", "<- wAB |"
      end if
    end if

    ! Let not have these other arrays hanging around a moment longer than we need them
    if (on_root) then
      deallocate (average_bandgap, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating average_bandgap in dos_utils: compute_bandgap')
      deallocate (optical_bandgap, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating optical_bandgap in dos_utils: compute_bandgap')
      deallocate (kpoint_bandgap, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating kpoint_bandgap in dos_utils: compute_bandgap')
      deallocate (thermal_vbm_multiplicity, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating thermal_vbm_multiplicity in dos_utils: compute_bandgap')
      deallocate (thermal_cbm_multiplicity, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocatingthermal_cbm_multiplicity in dos_utils: compute_bandgap')
      deallocate (optical_multiplicity, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating optical_multiplicity in dos_utils: compute_bandgap')
    end if

    if (on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) write (stdout, '(1x,a40,f11.3,a)') 'Time to calculate Bandgap ', time1 - time0, ' (sec)'
  end subroutine dos_utils_compute_bandgap
  !===============================================================================

  !===============================================================================
  function calc_band_energies(dos)
    !===============================================================================
    ! Function to return the band energy of a DOS by summing DOS*Energy up to
    ! Fermi energy
    !-------------------------------------------------------------------------------
    ! Arguments: dos (in) : Density of States array
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: none beyond the parent module variables
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: efermi must be set in od_electonic
    !                       E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries:
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !===============================================================================
    use od_electronic, only: nspins
    use od_parameters, only: dos_nbins
    use od_electronic, only: efermi
    implicit none

    real(kind=dp), intent(in) :: dos(1:dos_nbins, 1:nspins)
    real(kind=dp) :: gband, calc_band_energies
    integer :: is, idos

    gband = 0.0_dp
    do is = 1, nspins
      do idos = 1, dos_nbins
        if (E(idos) .le. efermi) then
          gband = gband + dos(idos, is)*E(idos)*delta_bins
        else
          exit
        end if
      end do
    end do

    calc_band_energies = gband
  end function calc_band_energies

  !===============================================================================
  function calc_efermi_from_intdos(INTDOS)
    !===============================================================================
    ! Function to calculate the Fermi energy from an intgrated DOS by serching for
    ! the bin which contains the correct number of electrons
    !-------------------------------------------------------------------------------
    ! Arguments: intdos (in) : Density of States array
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries:
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Modified from LinDOS
    !===============================================================================
    use od_electronic, only: num_electrons, nspins
    use od_parameters, only: dos_nbins
    use od_io, only: stdout
    use od_comms, only: on_root

    implicit none
    real(kind=dp), intent(in) :: INTDOS(1:dos_nbins, 1:nspins)
    real(kind=dp) :: calc_efermi_from_intdos
    real(kind=dp) :: efermi
    real(kind=dp) :: tolerance = 0.0001_dp ! Has the effect of forcing the efermi to
    ! be in the middle of the band gap, and not
    ! the leading edge.
    integer  :: idos, i, j

    do idos = 1, dos_nbins
      if (sum(INTDOS(idos, :)) .gt. (sum(num_electrons(:)) + (tolerance/2.0_dp))) exit
    end do

    efermi = E(idos)

    do i = idos, 1, -1
      if (sum(INTDOS(i, :)) .lt. (sum(INTDOS(idos - 1, :)) - (tolerance/2.0_dp))) exit
    end do

    calc_efermi_from_intdos = (efermi + E(i))/2.0_dp

    if (on_root) then
      do j = 1, nspins
        write (stdout, '(1x,a1,a20,i1,a20,f10.5,a4,1x,f10.5,3x,a8)') "|", "Spin Component : ", j,&
            &" occupation between ", INTDOS(i, j), "and", INTDOS(idos, j), "<- Occ |"
      end do
    end if

  end function calc_efermi_from_intdos

  !===============================================================================
  subroutine dos_utils_compute_band_energies
    !===============================================================================
    ! High-level subroutine to compute band energies of the DOS calculated.
    ! Calculates using the band_energies directly and compares with the
    ! function calc_band_energies which does the low level computation on the DOS.
    !-------------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E,dos_fixed,dos_adaptive,dos_linear
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !===============================================================================
    use od_parameters, only: adaptive, linear, fixed, iprint
    use od_electronic, only: electrons_per_state, efermi, nbands, nspins, band_energy
    use od_cell, only: kpoint_weight, num_kpoints_on_node
    use od_comms, only: comms_reduce, my_node_id, on_root
    use od_io, only: stdout, io_time

    implicit none
    real(kind=dp) :: eband
    real(kind=dp) :: time0, time1
    integer :: ik, is, ib

    time0 = io_time()
    if (on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Band Energy Analysis ---------------------------+'

      if (fixed) then
        write (stdout, '(1x,a1,a46,f12.4,a3,8x,a8)') "|", &
          " Band energy (Fixed broadening)  : ", calc_band_energies(dos_fixed), "eV", "<- BEF |"
      end if
      if (adaptive) then
        write (stdout, '(1x,a1,a46,f12.4,a3,8x,a8)') "|", &
          " Band energy (Adaptive broadening) : ", calc_band_energies(dos_adaptive), "eV", "<- BEA |"
      end if
      if (linear) then
        write (stdout, '(1x,a1,a46,f12.4,a3,8x,a8)') "|", &
          " Band energy (Linear broadening) : ", calc_band_energies(dos_linear), "eV", "<- BEL |"
      end if
    end if

    eband = 0.0_dp
    do ik = 1, num_kpoints_on_node(my_node_id)
      do is = 1, nspins
        do ib = 1, nbands
          if (band_energy(ib, is, ik) .le. efermi) eband = eband + band_energy(ib, is, ik)*electrons_per_state&
                                                          &*kpoint_weight(ik)
        end do
      end do
    end do
    call comms_reduce(eband, 1, 'SUM')
    if (on_root) then
      write (stdout, '(1x,a1,a46,f12.4,a3,8x,a8)') "|", " Band energy (From CASTEP) : ", eband, "eV", "<- BEC |"
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      time1 = io_time()
      if (iprint > 1) then
        write (stdout, '(1x,a59,f11.3,a8)') &
             '+ Time to calculate Band Energies                        &
          &      ', time1 - time0, ' (sec) +'
      end if
      write (stdout, '(1x,a78)')
    end if

  end subroutine dos_utils_compute_band_energies

  !===============================================================================
  subroutine setup_energy_scale
    !===============================================================================
    ! Sets up all broadening independent DOS concerns. That is basically the energy
    ! scale, E, and the width of its bins, delta_bins
    !-------------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: band_energy in od_electronic must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !===============================================================================
    use od_parameters, only: dos_nbins, dos_min_energy, dos_max_energy, dos_spacing, iprint
    use od_electronic, only: band_energy
    use od_io, only: io_error, stdout
    use od_comms, only: comms_reduce, comms_bcast, on_root

    implicit none

    real(kind=dp) :: min_band_energy, max_band_energy
    integer       :: idos, ierr

    if (allocated(E)) then
      deallocate (E, stat=ierr)
      if (ierr /= 0) call io_error("cannot deallocate E in dos_utils setup_energy_scale")
    end if

    ! If we do have dos_min_energy and dos_max_energy set, then we'd better
    ! use them. If not, let's set some sensible values.
    if (dos_min_energy == -huge(dos_min_energy)) then !Do it automatically
      min_band_energy = minval(band_energy) - 5.0_dp
    else
      min_band_energy = dos_min_energy
    end if
    call comms_reduce(min_band_energy, 1, 'MIN')
    call comms_bcast(min_band_energy, 1)

    if (dos_max_energy == huge(dos_max_energy)) then !Do it automatically
      max_band_energy = maxval(band_energy) + 5.0_dp
    else
      max_band_energy = dos_max_energy
    end if
    call comms_reduce(max_band_energy, 1, 'MAX')
    call comms_bcast(max_band_energy, 1)

    ! If dos_nbins is set, then we'd better use that
    if (dos_nbins < 0) then ! we'll have to work it out
      dos_nbins = abs(ceiling((max_band_energy - min_band_energy)/dos_spacing))
      ! Now modify the max_band_energy
      if (on_root .and. (iprint > 2)) then
        write (stdout, *)
        write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
        write (stdout, '(1x,a40,f11.3,13x,a14)') '| max_band_energy (before correction) : ', max_band_energy, "<-- DOS Grid |"
      end if
      max_band_energy = min_band_energy + dos_nbins*dos_spacing
    end if

    allocate (E(1:dos_nbins), stat=ierr)
    if (ierr /= 0) call io_error("cannot allocate E in dos_utils setup_energy_scale")

    delta_bins = (max_band_energy - min_band_energy)/real(dos_nbins - 1, dp)
    do idos = 1, dos_nbins
      E(idos) = min_band_energy + real(idos - 1, dp)*delta_bins
    end do

    if (on_root .and. (iprint > 2)) then
      write (stdout, '(1x,1a,a39,e11.5,13x,a14)') '|', 'dos_min_energy : ', dos_min_energy, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,e11.5,13x,a14)') '|', 'dos_max_energy : ', dos_max_energy, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,f11.3,13x,a14)') '|', 'min_band_energy : ', min_band_energy, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,f11.3,13x,a14)') '|', 'max_band_energy : ', max_band_energy, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,i11,13x,a14)') '|', 'dos_nbins : ', dos_nbins, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,f11.3,13x,a14)') '|', 'dos_spacing : ', dos_spacing, "<-- DOS Grid |"
      write (stdout, '(1x,1a,a39,f11.3,13x,a14)') '|', 'delta_bins : ', delta_bins, "<-- DOS Grid |"
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine setup_energy_scale

  !===============================================================================
  subroutine dos_utils_merge(dos, weighted_dos)
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
    use od_comms, only: comms_reduce
    use od_electronic, only: nspins
    use od_parameters, only: dos_nbins

    implicit none
    real(kind=dp), intent(inout), allocatable, optional :: weighted_dos(:, :, :) ! bins.spins, orbitals
    real(kind=dp), allocatable, intent(inout) :: dos(:, :)

    call comms_reduce(dos(1, 1), nspins*dos_nbins, "SUM")

    if (present(weighted_dos)) call comms_reduce(weighted_dos(1, 1, 1), mw%nspins*dos_nbins*mw%norbitals, "SUM")

!    if(.not.on_root) then
!       if(allocated(dos)) deallocate(dos,stat=ierr)
!       if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate dos")
!       if(present(weighted_dos))  then
!          if(allocated(weighted_dos)) deallocate(weighted_dos,stat=ierr)
!          if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate weighted_dos")
!       end if
!    endif
  end subroutine dos_utils_merge

  !===============================================================================
  subroutine allocate_dos(dos, intdos, w_dos)
    !===============================================================================
    ! Allocate the dos, intdos and w_dos if necessary.
    !-------------------------------------------------------------------------------
    ! Arguments: dos    (inout)       : The Density of States
    !            intdos (inout)       : The Integrated DOS
    !            w_dos  (inout) (opt) : Weighted DOS
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
    use od_electronic, only: nspins
    use od_parameters, only: dos_nbins
    use od_io, only: io_error

    implicit none

    real(kind=dp), allocatable, intent(inout)  :: dos(:, :)
    real(kind=dp), allocatable, intent(inout)  :: intdos(:, :)

    real(kind=dp), intent(out), optional, allocatable    :: w_dos(:, :, :) ! bins.spins, orbitals

    integer :: ierr

    allocate (dos(dos_nbins, nspins), stat=ierr)
    if (ierr /= 0) call io_error("error in allocating dos")
    dos = 0.0_dp

    allocate (intdos(dos_nbins, nspins), stat=ierr)
    if (ierr /= 0) call io_error("error in allocating intdos")
    intdos = 0.0_dp

    if (present(w_dos)) then
      allocate (w_dos(dos_nbins, mw%nspins, mw%norbitals), stat=ierr)
      if (ierr /= 0) call io_error("error in allocating weighted_dos")
      w_dos = 0.0_dp
    end if
  end subroutine allocate_dos

  subroutine dos_utils_deallocate
    use od_io, only: io_error
    implicit none
    integer :: ierr

    if (allocated(dos_adaptive)) then
      deallocate (dos_adaptive, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate dos_adaptive')
    end if
    if (allocated(dos_fixed)) then
      deallocate (dos_fixed, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate dos_fixed')
    end if
    if (allocated(dos_linear)) then
      deallocate (dos_linear, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate dos_linear')
    end if
    if (allocated(intdos_adaptive)) then
      deallocate (intdos_adaptive, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate intdos_adaptive')
    end if
    if (allocated(intdos_fixed)) then
      deallocate (intdos_fixed, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate intdos_fixed')
    end if
    if (allocated(intdos_linear)) then
      deallocate (intdos_linear, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate intdos_fixed')
    end if
    if (allocated(E)) then
      deallocate (E, stat=ierr)
      if (ierr /= 0) call io_error('Error: dos_utils_deallocate - failed to deallocate E')
    end if
  end subroutine dos_utils_deallocate

  !===============================================================================
  subroutine calculate_dos(dos_type, dos, intdos, matrix_weights, weighted_dos)
    !===============================================================================
    ! Once everything is set up this is the main workhorse of the module.
    ! It accumulates the DOS and WDOS be looping over spins, kpoints and bands.
    !-------------------------------------------------------------------------------
    ! Arguments: dos           (out)       : The Density of States
    !                          (in)        : one of "a", "f", "l", "q"
    !            intdos        (out)       : The Integrated DOS
    !            matrix_weights(in)  (opt) : The weightings, such as LCAO for the
    !                                        weighted dos.
    !            weighted_dos  (out) (opt) : Weighted DOS
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: We use local adaptive, fixed and linear variable so that if multiple
    ! are set, it won't always appear to do the linear scheme.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_constants, only: sqrt_two
    use od_algorithms, only: gaussian, algorithms_erf
    use od_cell, only: kpoint_grid_dim, kpoint_weight, num_kpoints_on_node, recip_lattice
    use od_electronic, only: band_gradient, electrons_per_state, nbands, nspins, band_energy
    use od_parameters, only: adaptive_smearing, fixed_smearing, hybrid_linear &
         &, finite_bin_correction, iprint, dos_nbins, numerical_intdos, hybrid_linear_grad_tol&
         &, linear_smearing
    use od_io, only: stdout, io_error
    use od_comms, only: my_node_id, on_root

    implicit none

    integer :: ik, is, ib, idos, iorb, i, ierr
    real(kind=dp) :: adaptive_smearing_temp, dos_temp, cuml, intdos_accum, width
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4)
    real(kind=dp) :: sub_cell_length(1:3)
    real(kind=dp) :: temp_debug

    character(len=1), intent(in)                    :: dos_type

    real(kind=dp), intent(out), allocatable, optional :: weighted_dos(:, :, :)
    real(kind=dp), intent(in), optional :: matrix_weights(:, :, :, :)

    real(kind=dp), intent(out), allocatable :: dos(:, :), intdos(:, :)
    real(kind=dp), allocatable :: dos_smear(:, :)

    logical :: linear, fixed, adaptive, force_adaptive

    linear = .false.
    fixed = .false.
    adaptive = .false.

    select case (dos_type)
    case ("l")
      linear = .true.
    case ("a")
      adaptive = .true.
    case ("f")
      fixed = .true.
    case default
      call io_error(" ERROR : unknown dos_type in calculate_dos ")
    end select

    width = 0.0_dp

    if (linear .or. adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:), dp)/2.0_dp
    if (adaptive .or. hybrid_linear) then
      do i = 1, 3
        sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
      end do
      adaptive_smearing_temp = adaptive_smearing*sum(sub_cell_length)/3.0_dp
    end if

    if (fixed) width = fixed_smearing

    if (calc_weighted_dos) then
      call allocate_dos(dos, intdos, w_dos=weighted_dos)
    else
      call allocate_dos(dos, intdos)
    end if
    if (iprint > 1 .and. on_root) then ! This is to contain the Calculating k-points block
      write (stdout, '(1x,a78)') '+------------------------------ Calculate DOS -------------------------------+'
    end if
    do ik = 1, num_kpoints_on_node(my_node_id)
      if (iprint > 1 .and. on_root) then
        if (mod(real(ik, dp), 10.0_dp) == 0.0_dp) write (stdout, '(1x,a1,a28,i4,a3,i4,a14,7x,a17)') "|",&
             &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id), " on this node", "<-- DOS |"
      End if
      do is = 1, nspins
        do ib = 1, nbands
          if (linear .or. adaptive) grad(:) = band_gradient(ib, :, ik, is)
          ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
          ! adaptive smearing (and take advantage of FBCS if required).
          force_adaptive = .false.

          if (linear .and. .not. force_adaptive) call doslin_sub_cell_corners(grad, step, band_energy(ib, is, ik), EV)
          if (adaptive .or. force_adaptive) width = sqrt(dot_product(grad, grad))*adaptive_smearing_temp
          ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
          ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
          if (finite_bin_correction .and. (width < delta_bins)) width = delta_bins

          intdos_accum = 0.0_dp

          do idos = 1, dos_nbins
            ! The linear method has a special way to calculate the integrated dos
            ! we have to take account for this here.
            if (linear .and. .not. force_adaptive) then
              dos_temp = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), E(idos), cuml)*electrons_per_state*kpoint_weight(ik)
              intdos(idos, is) = intdos(idos, is) + electrons_per_state*kpoint_weight(ik)*cuml
            else
              dos_temp = gaussian(band_energy(ib, is, ik), width, E(idos))*electrons_per_state*kpoint_weight(ik)
              if (numerical_intdos) then
                intdos_accum = intdos_accum + dos_temp
                intdos(idos, is) = intdos(idos, is) + intdos_accum
              else ! Do it (semi)-analytically
                intdos(idos, is) = intdos(idos, is) + 0.5_dp*(1.0_dp + algorithms_erf((E(idos) - &
                     & band_energy(ib, is, ik))/(sqrt_two*width)))*electrons_per_state*kpoint_weight(ik)
              end if
            end if

            dos(idos, is) = dos(idos, is) + dos_temp

            if (calc_weighted_dos) then
              if (ik .le. mw%nkpoints) then
                if (ib .le. mw%nbands) then
                  do iorb = 1, mw%norbitals
                    weighted_dos(idos, is, iorb) = weighted_dos(idos, is, iorb) + &
                         & dos_temp*matrix_weights(iorb, ib, ik, is)
                  end do
                end if
              end if
            end if
          end do
        end do
      end do
    end do

    if (linear .and. (linear_smearing .gt. 0.0_dp)) then
      if (iprint > 1 .and. on_root) then ! This is to contain the Calculating k-points block
        write (stdout, '(1x,a78)') '+-------------------------------- Smear DOS ---------------------------------+'
      end if
      ! Post smear the dos with a Guassian
      ! allocate a temporary array
      allocate (dos_smear(dos_nbins, nspins), stat=ierr)
      if (ierr /= 0) call io_error("error in allocating dos_smear")
      dos_smear = 0.0_dp
      ! loop over spins
      ! loop over temporary array
      ! loop over dos array
      do is = 1, nspins
        do idos = 1, dos_nbins
          if (iprint > 1 .and. on_root) then
            if (mod(real(idos, dp), 1000.0_dp) == 0.0_dp) write (stdout, '(1x,a1,a25,i14,a3,i10,a14,1x,a10)') "|",&
           &"Calculating bin ", idos, " of", dos_nbins, " on this node", "<--sDOS |"
          end if
          do i = 1, dos_nbins
            dos_smear(idos, is) = dos_smear(idos, is) + dos(i, is)*gaussian(E(idos), linear_smearing, E(i))*delta_bins
          end do
        end do
      end do
      ! copy array back
      dos = dos_smear
      ! deallocate array
      deallocate (dos_smear)
    end if

    if (iprint > 1 .and. on_root) then ! This is to contain the Calculating k-points block
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    if (.not. linear .and. numerical_intdos) intdos = intdos*delta_bins

  end subroutine calculate_dos

  !===============================================================================
  subroutine doslin_sub_cell_corners(grad, step, energy, EigenV)
    !===============================================================================
    ! A helper subroutine for calculated_dos, which is used for the linear
    ! broadening method. This routine extrapolates the energies at the corner of the
    ! sub cells by using the gradient at the centre of the cell
    !-------------------------------------------------------------------------------
    ! Arguments: grad   (in) : The Gradient of the band at the centre of a sub-cell
    !            step   (in) : The directions to the edge of the sub_cell
    !            energy (in) : The Band energy at the centre of the sub cell
    !            EigenV (out): The Energies at the corners of the sub-cell
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: None
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_algorithms, only: heap_sort
    use od_cell, only: recip_lattice

    implicit none
    real(kind=dp), intent(in)  :: grad(1:3), step(1:3)
    real(kind=dp), intent(out) :: EigenV(0:4)
    real(kind=dp), intent(in)  :: energy

    integer :: m, n, o, nn, i
    real(kind=dp) :: stepp(1:3), DE(1:8)

    nn = 0
    do m = -1, 1, 2
      do n = -1, 1, 2
        do o = -1, 1, 2
          nn = nn + 1
          stepp(1) = step(1)*m
          stepp(2) = step(2)*n
          stepp(3) = step(3)*o

          ! Reciprocal lattice in inverse Ang
          stepp = matmul(recip_lattice, stepp)

          ! DE in eV
          DE(nn) = dot_product(grad, stepp)

        end do
      end do
    end do

    ! Yes this is a hack to the sorter to work the right way around.
    DE = -DE
    call heap_sort(8, DE)
    DE = -DE

    ! WHY ARE WE STORING ALL THIS?
    !EV(0,ib,is,ik) = band_energy(ib,is,ik)
    !EV(1,ib,is,ik) = EV(0,ib,is,ik) + DE(5)
    !EV(2,ib,is,ik) = EV(0,ib,is,ik) + DE(6)
    !EV(3,ib,is,ik) = EV(0,ib,is,ik) + DE(7)
    !EV(4,ib,is,ik) = EV(0,ib,is,ik) + DE(8)
    EigenV(0) = Energy
    do i = 1, 4
      EigenV(i) = EigenV(0) + DE(4 + i)
    end do
  end subroutine doslin_sub_cell_corners

  !===============================================================================
  function doslin(e0, e1, e2, e3, e4, e, int)
    !===============================================================================
    ! Return the DoS contribution for a linear band portion and a cubic cell
    !-------------------------------------------------------------------------------
    ! Arguments: e0 (in) : Energy at centre of sub cell
    !   e1,e2,e3,e4 (in) : Energies of the four lowest corners of the sub cell
    !            e  (in) : Energy at which DOS is evaluated
    !          int  (out): Integrated DOS contribution for energy, E)
    ! (The function itself returns the DOS couribution for energy, E)
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: None
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : C J Pickard. From LinDOS. Extra Comments A J Morris Sept 2010
    !===============================================================================
    implicit none

    ! - E X T R A   C O M M E N T S   B Y   A J M
    real(kind=dp), intent(in)  :: e0, e1, e2, e3, e4
    real(kind=dp), intent(in)  :: e
    real(kind=dp), intent(out) :: int

    real(kind=dp) :: doslin
    real(kind=dp) :: alpha, e_use
    logical        :: shift

    ! ** Check the input arguments

    ! The inputs must be in ascending order. These are the four lowest corners of the
    ! cube (sub cell) around the k-point. E_0 is the energy of the k-point which would
    ! in a Gaussian smearing scheme, be the only info we have about this cube of recip
    ! space.
    if (.not. ((e4 .le. e3) .and. (e3 .le. e2) .and. (e2 .le. e1) .and. (e1 .le. e0))) then
      stop 'doslin: input arguments incorrect'
    end if

    ! ** Return if outside the range

    ! If the energy is below the smallest corner, then return. The CES doesn't cut this
    ! sub cell
    if (e .le. e4) then
      doslin = 0.0_dp ! TRUE
      int = 0.0_dp ! TRUE
      return
    end if

    ! Since the energy in the sub cell is linearly extrapolated, if the energy is twice as
    ! big as the energy difference between the smallest corner and the middle, then this
    ! energy is outdside the top of the cell, and whilst this cell doesn't contribute to the
    ! CES, it does to the occupation.
    if (e .ge. (2.0_dp*e0 - e4)) then
      doslin = 0.0_dp ! TRUE
      int = 1.0_dp ! TRUE
      return
    end if

    ! ** Special treatment if all vertices at the same energy
    ! If the CES is perfectly flat in the cell, then we don't want to be dividing by zero.
    ! the below just catches this case and forces alpha to be 1. This is fine as the else
    ! block below looks for the same problem and the answer comes out correctly to
    if ((e1 == e2) .and. (e1 == e3) .and. (e3 == e4)) then
      alpha = 1.0_dp
    else
      alpha = 1.0_dp/((e1 - e3)*(e1 - e4) + (e3 - e4)**2/3.0_dp - (e1 - e2)**2/3.0_dp + (e1 + e2 - e3 - e4)*(e0 - e1)) ! TRUE
    end if

    ! ** Flip if above e0

    ! If e is greater than the energy of the k-point, then we're going to subtract the
    ! contribution from a full cell, rather than add it to an empty one. The extrapolation
    ! is linear, so we can do this fine.
    if (e .gt. e0) then
      e_use = 2.0_dp*e0 - e ! TRUE
      shift = .true.
    else
      e_use = e ! TRUE
      shift = .false.
    end if

    ! ** The analytic constributions to the DOS and integrated DOS

    if (e_use .le. e4) then
      ! P O S S I B L E   S P E E D   U P
      ! If we ended up in here, something went wrong as this should already have been trapped.
      doslin = 0.0_dp ! TRUE
      int = 0.0_dp ! TRUE
    else if (e_use .lt. e3) then
      ! There isn't a problem with divide by zero here. Since if e3=e4 and e_use < e3
      ! we would have been caugh in the above if.
      doslin = (e_use - e4)**2/(e3 - e4)/2.0_dp ! TRUE
      int = (e_use - e4)**3/(e3 - e4)/6.0_dp ! TRUE
    else if (e_use .lt. e2) then
      doslin = (e_use - (e3 + e4)/2.0_dp) ! TRUE
      int = ((e3 - e4)**2/3.0_dp + (e_use - e3)*(e_use - e4))/2.0_dp ! TRUE
    else if (e_use .lt. e1) then
      doslin = (e1 + e2 - e3 - e4)/2.0_dp ! TRUE
      ! Ok, so the IF costs more than the maths. But this way also catches the
      ! divide by zero. Clever!
      if (e1 .ne. e2) doslin = doslin - (e1 - e_use)**2/(e1 - e2)/2.0_dp ! TRUE
      int = ((e2 - e4)*(e_use - e3) + (e1 - e3)*(e_use - e2) + (e3 - e4)**2/3.0_dp + &
             ((e1 - e_use)**3 - (e1 - e2)**3)/(e1 - e2)/3.0_dp)/2.0_dp ! TRUE
    else if (e_use .le. e0) then
      ! Check to see if the band is flat.
      if ((e1 + e2 - e3 - e4) .gt. 0.0_dp) then
        doslin = (e1 + e2 - e3 - e4)/2.0_dp ! TRUE
        int = ((e1 - e3)*(e1 - e4) + (e3 - e4)**2/3.0_dp - (e1 - e2)**2/3.0_dp + &
               (e1 + e2 - e3 - e4)*(e_use - e1))/2.0_dp ! TRUE
      else ! This can only happen if e1=e2=e3=e4,
        ! in this case we stop doing what we were doing and calculate the contirbution from the
        ! gradient simplistically. grad E = e0-e4.
        doslin = 1.0_dp/(e0 - e4)/2.0_dp ! TRUE
        int = (e_use - e4)/(e0 - e4)/2.0_dp ! SO YES, THIS DOES INTEGRATE FROM THE LINE ABOVE
      end if
    else
      write (*, *) e_use, e
      write (*, *) e0, e1, e2, e3, e4
      stop 'Got here, but not supposed to!'
    end if

    ! ** Normalise

    doslin = doslin*alpha

    if (shift) then
      int = 1.0_dp - int*alpha ! TRUE
    else
      int = int*alpha ! TRUE
    end if

    return
  end function doslin

  !===============================================================================
  subroutine dos_utils_calculate_at_e(energy, dos_at_e, matrix_weights, weighted_dos_at_e)
    !===============================================================================
    ! Main routine in dos module, drives the calculation of density of states for
    ! both task : dos and also if it is required elsewhere.
    !-------------------------------------------------------------------------------
    ! Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
    !            weighted_dos   (out)(opt) : Output DOS weigthed by matrix_weights
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
    ! intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
    ! efermi_linear, delta_bins, calc_weighted_dos
    !-------------------------------------------------------------------------------
    ! Modules Used: see below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: One of linear, adaptive or fixed must be .true.
    !-------------------------------------------------------------------------------
    ! Known Worries: (1) If more than one of linear, adaptive or fixed are set it
    ! uses the most complicated method.
    ! (2) It should be possible to pass optioinal arguments to sub programs as
    ! optional argumnets without checking whether they are there or not. g95 will
    ! allow this behaviour. gfotran will not.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !===============================================================================
    use od_io, only: stdout, io_time, io_error
    use od_comms, only: on_root, my_node_id
    use od_electronic, only: band_gradient, nspins, elec_read_band_gradient
    use od_parameters, only: linear, adaptive, fixed, quad, iprint
    use od_cell, only: nkpoints, num_kpoints_on_node
    implicit none

    !-------------------------------------------------------------------------------
    ! I N T E R N A L   V A R I A B L E S
    real(kind=dp) :: time0, time1
    real(kind=dp), intent(in), allocatable, optional  :: matrix_weights(:, :, :, :)
    real(kind=dp), intent(out), optional  :: weighted_dos_at_e(:, :) ! spins, orbitals RJN3Jun changed
    real(kind=dp), intent(in) :: energy
    real(kind=dp), intent(out) :: dos_at_e(1:3, nspins) ! fixed, adaptive, linear : spins
    !-------------------------------------------------------------------------------

    if (.not. (linear .or. adaptive .or. fixed .or. quad)) call io_error(" DOS: No Broadening Set")

    calc_weighted_dos = .false.
    if (present(matrix_weights)) calc_weighted_dos = .true.

    if (calc_weighted_dos) then
      mw%norbitals = size(matrix_weights, 1)
      mw%nbands = size(matrix_weights, 2)
      mw%nkpoints = size(matrix_weights, 3)
      mw%nspins = size(matrix_weights, 4)
    end if

    if (calc_weighted_dos) then
      if (mw%nspins .ne. nspins) call io_error("ERROR : DOS :  mw%nspins not equal to nspins.")
      if (mw%nkpoints .ne. num_kpoints_on_node(my_node_id)) &
           & call io_error("ERROR : DOS : mw%nkpoints not equal to nkpoints.")
    end if

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if (quad .or. linear .or. adaptive) then
      if (.not. allocated(band_gradient)) call elec_read_band_gradient
    end if
    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   D O S
    ! Now everything is set up, we can perform the dos accumulation in parellel

    time0 = io_time()

    if (fixed) then
      if (calc_weighted_dos .and. (.not. adaptive) .and. (.not. linear)) then
        call calculate_dos_at_e("f", energy, dos_at_e(1, :), matrix_weights=matrix_weights, &
             &weighted_dos_at_e=weighted_dos_at_e)
        call dos_utils_merge_at_e(dos_at_e(1, :), weighted_dos_at_e=weighted_dos_at_e)
      else
        call calculate_dos_at_e("f", energy, dos_at_e(1, :))
        call dos_utils_merge_at_e(dos_at_e(1, :))
      end if
    end if
    if (adaptive) then
      if (calc_weighted_dos .and. (.not. linear)) then
        call calculate_dos_at_e("a", energy, dos_at_e(2, :), matrix_weights=matrix_weights, &
             &weighted_dos_at_e=weighted_dos_at_e)
        call dos_utils_merge_at_e(dos_at_e(2, :), weighted_dos_at_e=weighted_dos_at_e)
      else
        call calculate_dos_at_e("a", energy, dos_at_e(2, :))
        call dos_utils_merge_at_e(dos_at_e(2, :))
      end if
    end if
    if (linear) then
      if (calc_weighted_dos) then
        call calculate_dos_at_e("l", energy, dos_at_e(3, :), matrix_weights=matrix_weights, &
             &weighted_dos_at_e=weighted_dos_at_e)
        call dos_utils_merge_at_e(dos_at_e(3, :), weighted_dos_at_e=weighted_dos_at_e)
      else
        call calculate_dos_at_e("l", energy, dos_at_e(3, :))
        call dos_utils_merge_at_e(dos_at_e(3, :))
      end if
    end if
    ! if(quad) then

    !
    ! endif

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to calculate DOS at certain energy                &
           &      ', time1 - time0, ' (sec) +'
    end if

    !-------------------------------------------------------------------------------

  end subroutine dos_utils_calculate_at_e

  !===============================================================================
  subroutine calculate_dos_at_e(dos_type, energy, dos_at_e, matrix_weights, weighted_dos_at_e)
    !===============================================================================
    ! Once everything is set up this is the main workhorse of the module.
    ! It accumulates the DOS and WDOS be looping over spins, kpoints and bands.
    !-------------------------------------------------------------------------------
    ! Arguments: dos           (out)       : The Density of States
    !            intdos        (out)       : The Integrated DOS
    !            matrix_weights(in)  (opt) : The weightings, such as LCAO for the
    !                                        weighted dos.
    !            weighted_dos  (out) (opt) : Weighted DOS
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: If linear, fixed and adaptive are all set. It produces the
    ! linear result, no matter what was intended. This would need to be modified to
    ! take an optinal argument, which would force only one of the above to be set
    ! within the subroutine. Since linear, fixed and adaptive are only non-mutually
    ! exculsive when debugging (such things as efermi) this isn't a priority
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_algorithms, only: gaussian
    use od_cell, only: kpoint_grid_dim, kpoint_weight, num_kpoints_on_node, &
        & recip_lattice
    use od_electronic, only: band_gradient, electrons_per_state, nbands, nspins, band_energy
    use od_parameters, only: adaptive_smearing, fixed_smearing&
         &, finite_bin_correction, iprint, hybrid_linear_grad_tol, hybrid_linear
    use od_io, only: stdout, io_error
    use od_comms, only: my_node_id, on_root

    implicit none

    integer :: ik, ib, is, iorb, i
    real(kind=dp) :: dos_temp, cuml, intdos_accum, width, adaptive_smearing_temp
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4), sub_cell_length(1:3)

    character(len=1), intent(in)                    :: dos_type   !! RJN3JUN changed
!    character(len=1)                   :: dos_type   !! RJN3Jun changed
    real(kind=dp), intent(out), optional :: weighted_dos_at_e(:, :)  !! RJN3Jun changed
    real(kind=dp), intent(in), optional :: matrix_weights(:, :, :, :)
    real(kind=dp), intent(in) :: energy
    real(kind=dp), intent(out) :: dos_at_e(nspins)

    logical :: linear, fixed, adaptive, have_weighted_dos, force_adaptive

    have_weighted_dos = .false.
    if (present(weighted_dos_at_e)) have_weighted_dos = .true.

    linear = .false.
    fixed = .false.
    adaptive = .false.

!    if(fixed)       dos_type="f"     !! RJN3Jun added
!    if(adaptive)  dos_type="a"       !! RJN3Jun added
!    if(linear)      dos_type="l"     !! RJN3Jun added

    select case (dos_type)
    case ("l")
      linear = .true.
    case ("a")
      adaptive = .true.
    case ("f")
      fixed = .true.
    case default
      call io_error(" ERROR : unknown dos_type in calculate_dos ")
    end select

    width = 0.0_dp ! Just in case
    dos_at_e = 0.0_dp

    if (linear .or. adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:), dp)/2.0_dp
    if (adaptive .or. hybrid_linear) then
      do i = 1, 3
        sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
      end do
      adaptive_smearing_temp = adaptive_smearing*sum(sub_cell_length)/3
    end if

    if (fixed) width = fixed_smearing

    if (have_weighted_dos) weighted_dos_at_e = 0.0_dp

    if ((iprint > 1) .and. on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') "+------------------------- Calculate DOS at energy --------------------------+"
    end if

    do ik = 1, num_kpoints_on_node(my_node_id)
      if (iprint > 1 .and. on_root) then
        if (mod(real(ik, dp), 10.0_dp) == 0.0_dp) write (stdout, '(1x,1a,a28,i4,a3,i4,a14,10x,a14)') '|', &
             &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id), " on this node", "<-- DOS at E |"
      end if

      do is = 1, nspins

        do ib = 1, nbands
          if (.not. fixed) then
            if (linear .or. adaptive) grad(:) = band_gradient(ib, :, ik, is)
            ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
            ! adaptive smearing (and take advantage of FBCS if required).
            force_adaptive = .false.
            if (hybrid_linear .and. (hybrid_linear_grad_tol > sqrt(dot_product(grad, grad)))) force_adaptive = .true.
            if (linear .and. .not. force_adaptive) call doslin_sub_cell_corners(grad, step, band_energy(ib, is, ik), EV)
            if (adaptive .or. force_adaptive) width = sqrt(dot_product(grad, grad))*adaptive_smearing_temp
          end if
          ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
          ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
          if (finite_bin_correction) then ! Force the compiler to do this the right way around
            if (width < delta_bins) then
              width = delta_bins
            end if
          end if

          intdos_accum = 0.0_dp

          ! The linear method has a special way to calculate the integrated dos
          ! we have to take account for this here.
          if (linear .and. .not. force_adaptive) then
            dos_temp = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), energy, cuml)*electrons_per_state*kpoint_weight(ik)

          else
            dos_temp = gaussian(band_energy(ib, is, ik), width, energy)*electrons_per_state*kpoint_weight(ik)

          end if

          dos_at_e(is) = dos_at_e(is) + dos_temp

          if (have_weighted_dos) then
            if (ik .le. mw%nkpoints) then
              if (ib .le. mw%nbands) then
                do iorb = 1, mw%norbitals
                  weighted_dos_at_e(is, iorb) = weighted_dos_at_e(is, iorb) + &
                       & dos_temp*matrix_weights(iorb, ib, ik, is)
                end do
              end if
            end if
          end if
        end do
      end do
    end do

    if ((iprint > 1) .and. on_root) write (stdout, '(1x,a78)') &
         & "+----------------------------------------------------------------------------+"

  end subroutine calculate_dos_at_e

  !===============================================================================
  subroutine dos_utils_merge_at_e(dos, weighted_dos_at_e)
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
    use od_comms, only: comms_reduce
    use od_electronic, only: nspins

    implicit none
    real(kind=dp), intent(inout), optional :: weighted_dos_at_e(:, :) ! bins.spins, orbitals RJN3Jun changed
    real(kind=dp), intent(inout) :: dos(nspins)

    call comms_reduce(dos(1), nspins, "SUM")

    if (present(weighted_dos_at_e)) call comms_reduce(weighted_dos_at_e(1, 1), mw%nspins*mw%norbitals, "SUM")

!    if(.not.on_root) then
!       if(present(weighted_dos_at_e))  then
!          if(allocated(weighted_dos_at_e)) deallocate(weighted_dos_at_e,stat=ierr)
!          if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate weighted_dos")
!       end if
!    endif
  end subroutine dos_utils_merge_at_e

end module od_dos_utils
