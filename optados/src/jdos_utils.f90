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
! MODULE od_jdos  OptaDOS - Joint Density of States
! This is the module that contains all if the JDOS routines. It is used through
! the global jdos_calculate subroutine
!-------------------------------------------------------------------------------
! NB. It should be possible to pass optioinal arguments to sub programs as
! optional argumnets without checking whether they are there or not. g95 will
! allow this behaviour. gfotran will not.
!===============================================================================
module od_jdos_utils
  use od_constants, only: dp

  implicit none

  !-------------------------------------------------------------------------------
  ! P U B L I C   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: jdos_adaptive(:, :)
  real(kind=dp), allocatable, public, save :: jdos_fixed(:, :)
  real(kind=dp), allocatable, public, save :: jdos_linear(:, :)

  integer, save :: jdos_nbins

  real(kind=dp), allocatable, public, save :: E(:)

  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S
  public :: jdos_utils_calculate
  public :: setup_energy_scale
  !-------------------------------------------------------------------------------

  real(kind=dp), save                   :: delta_bins ! Width of bins
  logical :: calc_weighted_jdos
  integer, allocatable, save :: vb_max(:)
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine jdos_utils_calculate(matrix_weights, weighted_jdos)
    !===============================================================================
    ! Main routine in dos module, drives the calculation of Density of states for
    ! both task : dos and also if it is required elsewhere.
    !===============================================================================
    use od_parameters, only: linear, fixed, adaptive, quad, iprint, dos_per_volume, photo, photo_slab_volume,&
                            &jdos_max_energy, jdos_spacing
    use od_electronic, only: elec_read_band_gradient, band_gradient, nspins, efermi_set
    use od_comms, only: on_root, comms_bcast
    use od_io, only: stdout, io_error, io_time, seedname
    use od_cell, only: cell_volume
    use od_dos_utils, only: dos_utils_set_efermi

    implicit none
    !integer :: ierr
    real(kind=dp) :: time0, time1

    real(kind=dp), intent(out), allocatable, optional    :: weighted_jdos(:, :, :)  !I've added this
    real(kind=dp), intent(in), optional  :: matrix_weights(:, :, :, :, :)               !I've added this

    integer :: N_geom, is, idos, wjdos_unit = 25
    logical :: print_weighted_jdos = .false.

    calc_weighted_jdos = .false.
    if (present(matrix_weights)) calc_weighted_jdos = .true.

    if (calc_weighted_jdos .eqv. .false.) then ! We are called just to provide jdos.
      if (allocated(E)) then
        if (on_root .and. iprint > 1) write (stdout, *) " Already calculated jdos, so returning..."
        return  ! The jdos has already been calculated previously so just return.
      end if
    end if

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if (quad .or. linear .or. adaptive) then
      if (.not. allocated(band_gradient)) call elec_read_band_gradient
    end if
    !-------------------------------------------------------------------------------

    if (.not. efermi_set) call dos_utils_set_efermi

    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   J D O S
    ! Now everything is set up, we can perform the dos accumulation in parellel
    time0 = io_time()

    call setup_energy_scale(E)

    if (fixed) then
      if (calc_weighted_jdos) then
        call calculate_jdos('f', jdos_fixed, matrix_weights, weighted_jdos)
        call jdos_utils_merge(jdos_fixed, weighted_jdos)
      else
        call calculate_jdos('f', jdos_fixed)
        call jdos_utils_merge(jdos_fixed)
      end if

    end if
    if (adaptive) then
      if (calc_weighted_jdos) then
        call calculate_jdos('a', jdos_adaptive, matrix_weights, weighted_jdos)
        call jdos_utils_merge(jdos_adaptive, weighted_jdos)
      else
        call calculate_jdos('a', jdos_adaptive)
        call jdos_utils_merge(jdos_adaptive)
      end if
    end if
    if (linear) then
      if (calc_weighted_jdos) then
        call calculate_jdos('l', jdos_linear, matrix_weights, weighted_jdos)
        call jdos_utils_merge(jdos_linear, weighted_jdos)
      else
        call calculate_jdos('l', jdos_linear)
        call jdos_utils_merge(jdos_linear)
      end if
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
      write (stdout, '(1x,a59,f11.3,a8)') &
           '+ Time to calculate Joint Density of States              &
           &      ', time1 - time0, ' (sec) +'
    end if
    !-------------------------------------------------------------------------------
    if (print_weighted_jdos) then
      if (on_root) then
        N_geom = size(matrix_weights, 5)
        open (unit=wjdos_unit, action='write', file=trim(seedname)//'_weighted_jdos.dat')
        write (wjdos_unit, '(1x,a28)') '############################'
        write (wjdos_unit, '(1x,a19,1x,a99)') '# Weighted JDOS for', seedname
        write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# maximum JDOS energy :', jdos_max_energy, '[eV]'
        write (wjdos_unit, '(1x,a23,1x,F10.4,1x,a4)') '# JDOS step size      :', jdos_spacing, '[eV]'
        write (wjdos_unit, '(1x,a28)') '############################'
        do is = 1, nspins
          write (wjdos_unit, *) 'Spin Channel :', is
          do idos = 1, jdos_nbins
            write (wjdos_unit, *) idos*jdos_spacing, ' , ', sum(weighted_jdos(idos, is, 1:N_geom))
          end do
        end do
        close (unit=wjdos_unit)
      end if
    end if

    if (dos_per_volume) then
      if (photo) then
        if (fixed) then
          jdos_fixed = jdos_fixed/photo_slab_volume
        end if
        if (adaptive) then
          jdos_adaptive = jdos_adaptive/photo_slab_volume
        end if
        if (linear) then
          jdos_linear = jdos_linear/photo_slab_volume
        end if
      else
        if (fixed) then
          jdos_fixed = jdos_fixed/cell_volume
        end if
        if (adaptive) then
          jdos_adaptive = jdos_adaptive/cell_volume
        end if
        if (linear) then
          jdos_linear = jdos_linear/cell_volume
        end if
      end if

      ! if(quad) then
      !    dos_quad=dos_quad/cell_volume
      !    intdos_quad=intdos_quad/cell_volume
      ! endif
    end if

  end subroutine jdos_utils_calculate

  !===============================================================================
  subroutine setup_energy_scale(E)
    !===============================================================================
    ! Sets up all broadening independent DOS concerns
    ! Calls the relevant dos calculator.
    !===============================================================================
    use od_dos_utils, only: dos_utils_calculate
    use od_parameters, only: jdos_max_energy, jdos_spacing, iprint
    use od_electronic, only: efermi, band_energy
    use od_comms, only: comms_reduce, comms_bcast, on_root
    use od_io, only: stdout, io_error

    implicit none

    integer       :: idos, ierr
    real(kind=dp) :: max_band_energy
    real(kind=dp), intent(out), allocatable, optional    :: E(:)

    if (jdos_max_energy < 0.0_dp) then ! we have to work it out ourselves
      max_band_energy = maxval(band_energy)
      call comms_reduce(max_band_energy, 1, 'MAX')
      call comms_bcast(max_band_energy, 1)
      jdos_max_energy = efermi - max_band_energy

      if (on_root .and. (iprint > 2)) then
        write (stdout, *)
        write (stdout, '(1x,a78)') &
             & '+----------------------------------------------------------------------------+'
        write (stdout, '(1x,a1,a38,f11.3,13x,a15)') &
             & '|', 'max_band_energy (before correction) : ',&
             & max_band_energy, "<-- JDOS Grid |"
      end if
    end if

    jdos_nbins = abs(ceiling(jdos_max_energy/jdos_spacing))
    jdos_max_energy = jdos_nbins*jdos_spacing

    allocate (E(1:jdos_nbins), stat=ierr)
    if (ierr /= 0) call io_error("Error: jdos_utils, setup_energy_scale: cannot allocate E")

    delta_bins = jdos_max_energy/real(jdos_nbins - 1, dp)
    do idos = 1, jdos_nbins
      E(idos) = real(idos - 1, dp)*delta_bins
    end do

    if (on_root .and. (iprint > 2)) then
      write (stdout, '(1x,a1,a38,f11.3,13x,a15)') '|', 'efermi : ', efermi, "<-- JDOS Grid |"
      write (stdout, '(1x,a1,a38,f11.3,13x,a15)') '|', 'jdos_max_energy : ', jdos_max_energy, "<-- JDOS Grid |"
      write (stdout, '(1x,a1,a38,i11,13x,a15)') '|', ' jdos_nbins : ', jdos_nbins, "<-- JDOS Grid |"
      write (stdout, '(1x,a1,a38,f11.3,13x,a15)') '|', 'jdos_spacing : ', jdos_spacing, "<-- JDOS Grid |"
      write (stdout, '(1x,a1,a38,f11.3,13x,a15)') '|', 'delta_bins : ', delta_bins, "<-- JDOS Grid |"
      write (stdout, '(1x,a78)') &
        '+----------------------------------------------------------------------------+'
    end if

  end subroutine setup_energy_scale

  !===============================================================================
  subroutine allocate_jdos(jdos)
    !===============================================================================
    !===============================================================================
    use od_electronic, only: nspins
    use od_io, only: io_error
    implicit none

    real(kind=dp), allocatable, intent(out)  :: jdos(:, :)

    integer :: ierr

    allocate (jdos(jdos_nbins, nspins), stat=ierr)
    if (ierr /= 0) call io_error("Error in allocating jdos (jdos_utils)")
    jdos = 0.0_dp

  end subroutine allocate_jdos

  !===============================================================================
  subroutine jdos_deallocate
    !===============================================================================
    !===============================================================================
    use od_io, only: io_error
    implicit none

    integer :: ierr

    if (allocated(jdos_adaptive)) then
      deallocate (jdos_adaptive, stat=ierr)
      if (ierr /= 0) call io_error('Error: jdos_deallocate - failed to deallocate jdos_adaptive')
    end if
    if (allocated(jdos_fixed)) then
      deallocate (jdos_fixed, stat=ierr)
      if (ierr /= 0) call io_error('Error: jdos_deallocate - failed to deallocate jdos_fixed')
    end if
    if (allocated(jdos_linear)) then
      deallocate (jdos_linear, stat=ierr)
      if (ierr /= 0) call io_error('Error: jdos_deallocate - failed to deallocate jdos_linear')
    end if
    if (allocated(E)) then
      deallocate (E, stat=ierr)
      if (ierr /= 0) call io_error('Error: jdos_deallocate - failed to deallocate E')
    end if

  end subroutine jdos_deallocate

  !===============================================================================
  subroutine calculate_jdos(jdos_type, jdos, matrix_weights, weighted_jdos)
    !===============================================================================

    !===============================================================================
    use od_comms, only: my_node_id, on_root
    use od_cell, only: num_kpoints_on_node, kpoint_grid_dim, kpoint_weight,&
         &recip_lattice
    use od_parameters, only: adaptive_smearing, fixed_smearing, iprint, &
         &finite_bin_correction, scissor_op, hybrid_linear_grad_tol, hybrid_linear, exclude_bands, num_exclude_bands
    use od_io, only: io_error, stdout
    use od_electronic, only: band_gradient, nbands, band_energy, nspins, electrons_per_state, &
         & efermi
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    implicit none

    integer :: ik, is, ib, idos, jb, i
    integer :: N2, N_geom, ierr
    real(kind=dp) :: dos_temp, cuml, width, adaptive_smearing_temp
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4), sub_cell_length(1:3)

    character(len=1), intent(in)                      :: jdos_type
    real(kind=dp), intent(inout), allocatable, optional :: weighted_jdos(:, :, :)
    real(kind=dp), intent(in), optional                :: matrix_weights(:, :, :, :, :)

    real(kind=dp), intent(out), allocatable :: jdos(:, :)

    logical :: linear, fixed, adaptive, force_adaptive

    linear = .false.
    fixed = .false.
    adaptive = .false.

    select case (jdos_type)
    case ("l")
      linear = .true.
    case ("a")
      adaptive = .true.
    case ("f")
      fixed = .true.
    case default
      call io_error(" ERROR : unknown jdos_type in calculate_jdos ")
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

    call allocate_jdos(jdos)
    if (calc_weighted_jdos) then
      N_geom = size(matrix_weights, 5)
      allocate (weighted_jdos(jdos_nbins, nspins, N_geom), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_jdos - failed to allocate weighted_jdos')
      weighted_jdos = 0.0_dp
    end if

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------------ Calculate JDOS ------------------------------+'
    end if

    do ik = 1, num_kpoints_on_node(my_node_id)
      if (iprint > 1 .and. on_root) then
        if (mod(real(ik, dp), 10.0_dp) == 0.0_dp) write (stdout, '(1x,a1,a38,i4,a3,i4,1x,a14,3x,a10)') ',', &
             &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id), 'on this node.', "<-- JDOS |"
      end if
      do is = 1, nspins
        occ_states: do ib = 1, nbands
          if (num_exclude_bands > 0) then
            if (any(exclude_bands == ib)) cycle
          end if
          if (band_energy(ib, is, ik) .ge. efermi) cycle occ_states
          unocc_states: do jb = 1, nbands
            if (band_energy(jb, is, ik) .lt. efermi) cycle unocc_states
            if (linear .or. adaptive) grad(:) = band_gradient(jb, :, ik, is) - band_gradient(ib, :, ik, is)
            ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
            ! adaptive smearing (and take advantage of FBCS if required).
            force_adaptive = .false.
            if (.not. fixed) then
              if (hybrid_linear .and. (hybrid_linear_grad_tol > sqrt(dot_product(grad, grad)))) force_adaptive = .true.
              if (linear .and. .not. force_adaptive) call doslin_sub_cell_corners(grad, step, band_energy(jb, is, ik) -&
  &band_energy(ib, is, ik) + scissor_op, EV)
              if (adaptive .or. force_adaptive) width = sqrt(dot_product(grad, grad))*adaptive_smearing_temp
            end if

            ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
            ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
            if (finite_bin_correction .and. (width < delta_bins)) width = delta_bins

            do idos = 1, jdos_nbins
              ! The linear method has a special way to calculate the integrated dos
              ! we have to take account for this here.
              if (linear .and. .not. force_adaptive) then
                dos_temp = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), E(idos), cuml)
              else
                dos_temp=gaussian(band_energy(jb,is,ik)-band_energy(ib,is,ik)+scissor_op,width,E(idos))!&
              end if

              jdos(idos, is) = jdos(idos, is) + dos_temp*electrons_per_state*kpoint_weight(ik)

              ! this will become a loop over final index (polarisation)
              ! Also need to remove kpoints weights.
              if (calc_weighted_jdos) then
                do N2 = 1, N_geom
                  weighted_jdos(idos, is, N2) = weighted_jdos(idos, is, N2) + dos_temp*matrix_weights(ib, jb, ik, is, N2)&
                                               &*electrons_per_state*kpoint_weight(ik)
                end do
              end if

            end do
          end do unocc_states
        end do occ_states
      end do
    end do

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

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
    use od_comms, only: comms_reduce
    use od_electronic, only: nspins

    implicit none
    real(kind=dp), intent(inout), allocatable, optional :: weighted_jdos(:, :, :) ! bins.spins, orbitals
    real(kind=dp), allocatable, intent(inout) :: jdos(:, :)

    integer :: N_geom
    if (present(weighted_jdos)) N_geom = size(weighted_jdos, 3)

    call comms_reduce(jdos(1, 1), nspins*jdos_nbins, "SUM")

    if (present(weighted_jdos)) call comms_reduce(weighted_jdos(1, 1, 1), nspins*jdos_nbins*N_geom, "SUM")

!    if(.not.on_root) then
!       if(allocated(jdos)) deallocate(jdos,stat=ierr)
!       if (ierr/=0) call io_error (" ERROR : jdos : merge_jdos : cannot deallocate dos")
!       if(present(weighted_jdos))  then
!          if(allocated(weighted_jdos)) deallocate(weighted_jdos,stat=ierr)
!          if (ierr/=0) call io_error (" ERROR : jdos : merge_jdos : cannot deallocate weighted_dos")
!       end if
!    endif
  end subroutine jdos_utils_merge

end module od_jdos_utils
