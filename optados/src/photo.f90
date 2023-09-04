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
module od_photo

  use od_constants, only: dp

  implicit none
  private
  public :: photo_calculate

  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: pdos_weights_atoms
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: projected_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :, :, :) :: foptical_matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :) :: weighted_jdos
  real(kind=dp), allocatable, public, dimension(:, :) :: absorp_layer
  real(kind=dp), allocatable, dimension(:, :, :) :: pdos_weights_k_band
  real(kind=dp), allocatable, dimension(:, :, :) :: imfp_val
  real(kind=dp), allocatable, dimension(:, :, :, :) :: electron_esc
  real(kind=dp), dimension(:, :), allocatable :: I_layer
  real(kind=dp), allocatable, public, save :: E(:)
  real(kind=dp), allocatable, dimension(:, :) :: reflect_photo
  real(kind=dp), allocatable, dimension(:, :) :: absorp_photo

  real(kind=dp), allocatable, dimension(:, :) :: refract
  real(kind=dp), allocatable, dimension(:)  :: reflect
  real(kind=dp), allocatable, dimension(:) :: absorp

  real(kind=dp), dimension(:), allocatable :: thickness_atom
  real(kind=dp), dimension(:), allocatable :: thickness_layer
  real(kind=dp), dimension(:), allocatable :: atom_imfp
  real(kind=dp), dimension(:, :), allocatable :: new_atoms_coordinates
  real(kind=dp), allocatable, dimension(:, :, :) :: phi_arpes
  real(kind=dp), allocatable, dimension(:, :, :) :: theta_arpes
  real(kind=dp), allocatable, dimension(:, :, :) :: E_kinetic
  real(kind=dp), allocatable, dimension(:, :, :) :: E_transverse
  real(kind=dp), allocatable, dimension(:, :, :) :: bulk_prob
  real(kind=dp), allocatable, dimension(:) :: t_energy
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: weighted_temp
  integer :: max_energy
  real(kind=dp), allocatable, dimension(:, :, :, :) :: qe_osm
  real(kind=dp), allocatable, dimension(:, :, :, :, :) :: qe_tsm
  real(kind=dp) :: mean_te
  real(kind=dp) :: total_qe
  real(kind=dp), allocatable, dimension(:) :: layer_qe
  integer, dimension(:), allocatable :: atom_order
  integer, dimension(:), allocatable :: atoms_per_layer
  real(kind=dp) :: work_function_eff
  real(kind=dp) :: evacuum
  real(kind=dp) :: evacuum_eff
  real(kind=dp) :: total_field_emission
  real(kind=dp), allocatable, dimension(:, :, :) :: field_emission
  integer, allocatable, dimension(:) :: layer
  integer :: N_geom
  integer :: max_atoms
  integer :: max_layer
  real(kind=dp) :: q_weight

  ! Added by Felix Mildner, 12/2022

  integer, allocatable, dimension(:)  :: index_energy
  integer                             :: number_energies, current_energy_index, current_index
  real(kind=dp)                       :: temp_photon_energy, time_a, time_b
  integer, allocatable, dimension(:, :):: min_index_unocc
contains

  subroutine photo_calculate
    !
    !  Program to calculate the photoemission
    !

    use od_electronic, only: elec_dealloc_optical, elec_pdos_read, elec_read_optical_mat, &
    & efermi, efermi_set, elec_read_foptical_mat, elec_dealloc_pdos
    use od_jdos_utils, only: jdos_utils_calculate, setup_energy_scale
    use od_comms, only: on_root
    use od_parameters, only: photo_work_function, photo_model, photo_elec_field, write_photo_output, photo_photon_sweep,&
                            &photo_photon_min, jdos_spacing, photo_photon_energy, iprint, devel_flag
    use od_dos_utils, only: dos_utils_set_efermi, dos_utils_calculate_at_e, dos_utils_deallocate
    use od_io, only: stdout, io_error, io_time
    use od_pdos, only: pdos_calculate

    implicit none

    integer :: i

    if (on_root) then
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                             Photoemission Calculation                      +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    if (index(devel_flag, 'geom_analysis') > 0) then
      write (stdout, '(1x,a78)') '+           Only performing the analysis of the supplied geometry!           +'
      call calc_layers
      return
    end if

    if (.not. efermi_set) then
      call dos_utils_set_efermi
      call dos_utils_deallocate
    end if

    !Identify layers

    call calc_layers
    call calc_band_info

    call elec_read_optical_mat
    !THIS PART COMES FROM THE PDOS MODULE
    ! read in the pdos weights
    call elec_pdos_read
    call make_pdos_weights_atoms
    call elec_dealloc_pdos

    ! Calculate the optical properties of the slab
    call calc_photo_optics

    call calc_absorp_layer

    !Electric field and field emission
    if (photo_elec_field .gt. 0.0_dp) then
      call effect_wf
      call calc_field_emission
    else
      evacuum_eff = efermi + photo_work_function
      work_function_eff = photo_work_function
    end if

    if (photo_photon_sweep) then
      do i = 1, number_energies
        time_a = io_time()
        temp_photon_energy = photo_photon_min + (i - 1)*jdos_spacing
        if (on_root) write (stdout, '(1x,a50,f8.4,a20)') '+--------------- Starting Photoemission Sweep with', temp_photon_energy,&
                &' eV ---------------+'
        current_index = i
        current_energy_index = index_energy(i)
        !Calculate the photoemission angles theta/phi and transverse energy
        call calc_angle

        !Calculate the electron escape length
        call calc_electron_esc

        call bulk_emission

        !Calculate the QE
        if (index(photo_model, '3step') > 0) then !Three-step-model
          call calc_three_step_model
        elseif (index(photo_model, '1step') > 0) then !One-step-model
          if (.not. allocated(foptical_matrix_weights)) call elec_read_foptical_mat !Read the one-step matrix elements
          call make_foptical_weights !Calculate the one-step optical matrix
          call calc_one_step_model !Calculate QE
        end if

        !Weight the contribution of each electron
        !to the transverse energy spread according to their QE
        call weighted_mean_te

        call write_qe_data
        ! Only call the binding energy gaussian broadening and file printing if necessary
        if (.not. index(write_photo_output, 'off') > 0) then
          !Broaden ouputs using a gaussian function
          if (index(write_photo_output, 'e_bind') > 0) call binding_energy_spread
          !Write either a binding energy output with after Gaussian broadening or the reduced QE tensor
          call write_qe_output_files
        end if
        time_b = io_time()
        if (on_root .and. iprint > 1) then
          write (stdout, '(1x,a44,15x,f11.3,a8)') '+ Time to calculate Photoemission sweep step', time_b - time_a, ' (sec) +'
        end if
      end do
    else
      temp_photon_energy = photo_photon_energy
      current_index = 1
      current_energy_index = index_energy(1)
      !Calculate the photoemission angles theta/phi and transverse energy
      call calc_angle

      !Calculate the electron escape length
      call calc_electron_esc

      call bulk_emission

      !Calculate the QE
      if (index(photo_model, '3step') > 0) then !Three-step-model
        call calc_three_step_model
      elseif (index(photo_model, '1step') > 0) then !One-step-model
        if (.not. allocated(foptical_matrix_weights)) call elec_read_foptical_mat !Read the one-step matrix elements
        call make_foptical_weights !Calculate the one-step optical matrix
        call calc_one_step_model !Calculate QE
      end if

      !Weight the contribution of each electron
      !to the transverse energy spread according to their QE
      call weighted_mean_te
      call write_qe_data

      ! Only call the binding energy gaussian broadening and file printing if necessary
      if (.not. index(write_photo_output, 'off') > 0) then
        !Broaden ouputs using a gaussian function
        if (index(write_photo_output, 'e_bind') > 0) call binding_energy_spread
        !Write either a binding energy output with after Gaussian broadening
        call write_qe_output_files
      end if

    end if
    ! Deallocate the rest that was needed for the photoemission calcs
    call photo_deallocate

    if (on_root) write (stdout, '(1x,a78)') '| End of Photoemission Calculation                                           |'

  end subroutine photo_calculate

  !***************************************************************
  subroutine calc_layers
    !***************************************************************
    !This subroutine identifies the layer of each atom
    use od_constants, only: dp
    use od_cell, only: num_atoms, atoms_pos_cart_photo, atoms_label_tmp
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    use od_parameters, only: devel_flag, photo_max_layer, photo_layer_choice, photo_imfp_const
    implicit none
    integer :: atom_1, atom_2, i, atom_index, temp, first, ierr, atom, ic

    allocate (atom_order(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_layers - allocation of atom_order failed')

    allocate (layer(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_layers - allocation of layer failed')
    do i = 1, num_atoms
      atom_order(i) = i
    end do

    !SORTING ALGORITHM
    do atom_1 = 1, num_atoms - 1
      first = atom_order(atom_1)
      do atom_2 = atom_1 + 1, num_atoms
        atom_index = atom_1
        if (atoms_pos_cart_photo(3, atom_order(atom_2)) .gt. atoms_pos_cart_photo(3, first)) then
          first = atom_order(atom_2)
          atom_index = atom_2
        end if

        if (atom_index /= atom_1) then
          temp = atom_order(atom_1)
          atom_order(atom_1) = atom_order(atom_index)
          atom_order(atom_index) = temp
        end if
      end do
    end do

    ! DEFINE THE LAYER FOR EACH ATOM
    ! Assume that a new layer starts if the atom type changes or
    ! the atom is more than 0.5 Angstrom lower than the current layer
    i = 1
    layer(1) = 1
    do atom = 2, num_atoms
      if ((trim(atoms_label_tmp(atom_order(atom))) .ne. trim(atoms_label_tmp(atom_order(atom - 1)))) .or. &
          (abs(atoms_pos_cart_photo(3, atom_order(atom)) - atoms_pos_cart_photo(3, atom_order(atom - 1))) .gt. 0.50)) then
        i = i + 1
      end if
      layer(atom) = i
    end do

    do atom = 1, num_atoms
      ic = ichar(atoms_label_tmp(atom_order(atom)) (1:1))
      if ((ic .ge. ichar('a')) .and. (ic .le. ichar('z'))) &
        atoms_label_tmp(atom_order(atom)) (1:1) = char(ic + ichar('Z') - ichar('z'))
    end do
    if (on_root) then
      write (stdout, '(1x,a78)') '+------------------------------- Atomic Order  ------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Layer   |         Atom Z-Coordinate (Ang)          |'

      do atom = 1, num_atoms
        write (stdout, '(1x,a3,a2,8x,i3,11x,i3,18x,F12.7,a18)') "|  ", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
          layer(atom), &
          atoms_pos_cart_photo(3, atom_order(atom)), "|"
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    !CALCULATE THE MAX LAYER (HALF SLAB)
    max_layer = ((layer(num_atoms) + 1)/2)
    if (index(photo_layer_choice, 'user') > 0) then
      max_layer = photo_max_layer
    end if

    !CALCULATE THE MAX ATOM (HALF SLAB)
    max_atoms = 0
    do atom = 1, num_atoms
      if (layer(atom) .le. max_layer) then
        max_atoms = max_atoms + 1
      end if
    end do

    if (on_root) then
      write (stdout, 226) '|  Max number of atoms:', max_atoms, '   Max  number of layers:', max_layer, '   |'
226   format(1x, a23, I12, 1x, a25, 1x, I12, a4)

      if (index(devel_flag, 'layer_user') > 0) then
        write (stdout, '(1x,a78)') '|     *** ATTENTION *** : The max_layer value was supplied by the user!      |'
      end if

      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

    !CALCULATE HOW MANY ATOMS PER LAYER
    allocate (atoms_per_layer(max_layer), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_layers - allocation of atoms_per_layer failed')
    atoms_per_layer = 1
    do atom = 2, max_atoms
      if (layer(atom) .eq. layer(atom - 1)) then
        atoms_per_layer(layer(atom)) = atoms_per_layer(layer(atom)) + 1
      end if
    end do
    ! do i = 1, max_layer
    !   write (stdout, *) 'Layer: ', i, atoms_per_layer(i), ' |'
    ! end do

    !TEST IF THE SUPPLIED IMFP LIST IS LONG ENOUGH
    if (allocated(photo_imfp_const) .and. size(photo_imfp_const, 1) .gt. 1 .and. size(photo_imfp_const, 1) .lt. max_layer) then
      call io_error('The supplied list of layer dependent imfp values is less than the calculated max_layer. Check input!')
    end if

  end subroutine calc_layers

  subroutine calc_band_info
    !===============================================================================
    ! This subroutine determines useful indices of band energies for later use in
    ! the QE and MTE calculation to reduce loop times.
    ! This relies on an IMPORTANT assumption: the bands file is ordered by energy
    ! and not by band number (e.g. after being processed by bands2orbitals)
    ! Felix Mildner, 28th March 2023
    !===============================================================================

    use od_electronic, only: efermi, band_energy, nbands, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: iprint
    use od_io, only: stdout, io_time, io_error

    integer         :: N, N_spin, n_eigen, ierr
    real(kind=dp)   :: time0, time1

    time0 = io_time()

    allocate (min_index_unocc(nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_band_info - allocation of min_index_unocc failed')

    do N = 1, num_kpoints_on_node(my_node_id)  ! Loop over kpoints
      do N_spin = 1, nspins                           ! Loop over spins
        do n_eigen = 2, nbands                        ! Loop over bands
          ! TODO: Test if this is the behaviour we want and or if we have to change the condition
          if (band_energy(n_eigen - 1, N_spin, N) .gt. band_energy(n_eigen, N_spin, N)) then
            call io_error('Error: the band energies in the .bands file used are NOT ORDERED CORRECTLY (i.e. by increasing energy) &
            & which will give WRONG RESULTS!')
          end if
        end do
      end do
    end do

    do N = 1, num_kpoints_on_node(my_node_id)  ! Loop over kpoints
      do N_spin = 1, nspins                           ! Loop over spins
        do n_eigen = 1, nbands                        ! Loop over bands
          ! TODO: Test if this is the behaviour we want and or if we have to change the condition
          if (band_energy(n_eigen, N_spin, N) .gt. efermi) then
            min_index_unocc(N_spin, N) = n_eigen
            exit
          end if
        end do
      end do
    end do

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a36,23x,f11.3,a8)') '+ Time to calculate Band Energy Info', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_band_info

  !***************************************************************
  subroutine make_pdos_weights_atoms
    !***************************************************************
    !This subroutine is equivalent to pdos_merge of pdos.F90, but only for atoms
    use od_electronic, only: pdos_orbital, pdos_weights, pdos_mwab, nspins
    use od_cell, only: num_kpoints_on_node, num_atoms
    use od_comms, only: my_node_id, on_root
    use od_io, only: io_error, stdout
    use od_parameters, only: devel_flag
    implicit none
    integer :: N, N_spin, n_eigen, np, ierr, atom, i, i_max

    allocate (pdos_weights_atoms(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id), num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')

    allocate (pdos_weights_k_band(pdos_mwab%nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_k_band failed')

    pdos_weights_atoms = 0.0_dp
    pdos_weights_k_band = 0.0_dp

    do N = 1, num_kpoints_on_node(my_node_id)
      do N_spin = 1, nspins
        do n_eigen = 1, pdos_mwab%nbands
          i = 1
          do np = 1, pdos_mwab%norbitals
            if (np .gt. 1) then
              if (pdos_orbital%rank_in_species(np) .ne. pdos_orbital%rank_in_species(np - 1)) then
                i = i + 1
              end if
            end if
            pdos_weights_atoms(n_eigen, N_spin, N, i) = &
              pdos_weights_atoms(n_eigen, N_spin, N, i) + &
              pdos_weights(np, n_eigen, N, N_spin)
          end do
        end do
      end do
    end do
    i_max = i
    do atom = 1, max_atoms
      do N = 1, num_kpoints_on_node(my_node_id)
        do N_spin = 1, nspins
          do n_eigen = 1, pdos_mwab%nbands
            if (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom)) .lt. 0.0_dp) then
              pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom)) = 0.0_dp
            end if
            pdos_weights_k_band(n_eigen, N_spin, N) = pdos_weights_k_band(n_eigen, N_spin, N) + &
                                                      pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom))
          end do
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------ Printing pDOS_weights_atoms -----------------------+'
      write (stdout, 125) shape(pdos_weights_atoms)
      write (stdout, 125) i_max, pdos_mwab%nbands, num_kpoints_on_node(my_node_id), nspins
125   format(4(1x, I4))
      write (stdout, '(9999(es15.8))') ((((pdos_weights_atoms(n_eigen, N_spin, N, i), N_spin=1, nspins) &
                                          , n_eigen=1, pdos_mwab%nbands), N=1, num_kpoints_on_node(my_node_id)), i=1, i_max)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
      write (stdout, '(1x,a78)') '+----------------------- Printing pDOS_weights_k_band -----------------------+'
      write (stdout, 124) shape(pdos_weights_k_band)
      write (stdout, 124) pdos_mwab%nbands, num_kpoints_on_node(my_node_id), nspins
124   format(3(1x, I4))
      write (stdout, '(9999(es15.8))') (((pdos_weights_k_band(n_eigen, N_spin, N), &
                                          N=1, num_kpoints_on_node(my_node_id)), N_spin=1, nspins), n_eigen=1, pdos_mwab%nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if
  end subroutine make_pdos_weights_atoms

  !***************************************************************
  subroutine calc_photo_optics
    !***************************************************************

    use od_optics, only: make_weights, calc_epsilon_2, calc_epsilon_1, calc_refract, calc_absorp, calc_reflect, &
    & epsilon, refract, absorp, reflect, intra
    use od_io, only: stdout, io_error, io_time
    use od_electronic, only: elec_read_optical_mat, nbands, nspins, efermi, elec_dealloc_optical, elec_read_band_gradient,&
    & nbands, nspins, band_energy
    use od_cell, only: num_kpoints_on_node, num_kpoints_on_node
    use od_jdos_utils, only: jdos_utils_calculate, jdos_nbins, setup_energy_scale, jdos_deallocate
    use od_comms, only: comms_bcast, on_root, my_node_id
    use od_parameters, only: optics_intraband, jdos_spacing, photo_model, photo_photon_energy, photo_photon_sweep, &
      photo_photon_min, photo_photon_max, devel_flag, iprint
    use od_dos_utils, only: dos_utils_calculate_at_e
    use od_constants, only: epsilon_0, e_charge

    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: dos_matrix_weights
    real(kind=dp), allocatable, dimension(:, :) :: weighted_dos_at_e
    real(kind=dp), allocatable, dimension(:, :) :: dos_at_e

    integer :: N, N2, N_spin, n_eigen, n_eigen2, atom, ierr, energy
    integer :: jdos_bin, i, s
    real(kind=dp)    :: num_energies, temp, time0, time1

    time0 = io_time()

    if (photo_photon_sweep) then
      num_energies = (photo_photon_max - photo_photon_min)/jdos_spacing
      number_energies = int(num_energies) + 1
      if (photo_photon_max - photo_photon_min .eq. 0.0_dp) then
        number_energies = 1
      else if (mod(num_energies, 1.0_dp) > 1.0E-10_dp) then
        number_energies = number_energies + 1
        if (abs(mod(num_energies, 1.0_dp) - 1) > 1.0E-10_dp) &
          call io_error('Error: calc_photo_optics - the supplied photon sweep min/max values do not give integer # of photon steps')
      end if
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of index_energy failed')
      do i = 1, number_energies
        temp = (i - 1)*jdos_spacing + photo_photon_min
        ! Account for E = 0.0
        index_energy(i) = int(temp/jdos_spacing) + 1
      end do
    else
      number_energies = 1
      allocate (index_energy(number_energies), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of index_energy failed')
      ! Account for E = 0.0
      index_energy(number_energies) = int(photo_photon_energy/jdos_spacing) + 1
    end if

    allocate (absorp_photo(max_atoms, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of absorp_photo failed')

    allocate (reflect_photo(max_atoms, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of reflect_photo failed')

    call make_weights(matrix_weights)
    N_geom = size(matrix_weights, 5)
    call elec_dealloc_optical

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+-------------------------- Printing Matrix Weights -------------------------+'
      write (stdout, 126) shape(matrix_weights)
      write (stdout, 126) nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
      write (stdout, '(9999(es15.8))') (((((matrix_weights(n_eigen, n_eigen2, N, N_spin, N2), N2=1, N_geom), N_spin=1, nspins) &
                                          , N=1, num_kpoints_on_node(my_node_id)), n_eigen2=1, nbands), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    allocate (projected_matrix_weights(nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of projected_matrix_weights failed')

    do atom = 1, max_atoms                           ! Loop over atoms
      !
      if (iprint > 1 .and. on_root) then
        write (stdout, 145) '+------------------------ Starting Atom # ', atom, ' of ', max_atoms, ' ------------------------+'
145     format(1x, a42, I3, a4, I3, a26)
      end if
      ! (Re-)Setting the weights for new atom
      projected_matrix_weights = 0.0_dp

      do N2 = 1, N_geom
        do N = 1, num_kpoints_on_node(my_node_id)    ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands               ! Loop over state 1
              do n_eigen2 = n_eigen, nbands    ! Loop over state 2
                if (band_energy(n_eigen, N_spin, N) > efermi .and. n_eigen /= n_eigen2) cycle
                if (band_energy(n_eigen2, N_spin, N) < efermi .and. n_eigen /= n_eigen2) cycle
                projected_matrix_weights(n_eigen, n_eigen2, N, N_spin, N2) = &
                  matrix_weights(n_eigen, n_eigen2, N, N_spin, N2)* &
                  (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom))/pdos_weights_k_band(n_eigen, N_spin, N))
              end do                        ! Loop over state 2
            end do                            ! Loop over state 1
          end do                                ! Loop over spins
        end do                                    ! Loop over kpoints
      end do

      if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
        write (stdout, '(1x,a37,I3,a38)') '+-------------------------------Atom-', atom, '-------------------------------------+'
        write (stdout, '(1x,a78)') '+--------------------- Printing Projected Matrix Weights --------------------+'
        write (stdout, 126) shape(projected_matrix_weights)
        write (stdout, 126) nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
126     format(5(1x, I4))
        write (stdout, '(9999(es15.8))') (((((projected_matrix_weights(n_eigen, n_eigen2, N, N_spin, N2), N2=1, N_geom), &
                                             N_spin=1, nspins), N=1, num_kpoints_on_node(my_node_id)), &
                                           n_eigen2=1, nbands), n_eigen=1, nbands)
        write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
      end if

      ! Send matrix element to jDOS routine and get weighted jDOS back
      call jdos_utils_calculate(projected_matrix_weights, weighted_jdos)

      if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
        write (stdout, '(1x,a78)') '+------------------------ Printing Weighted Joint-DOS -----------------------+'
        write (stdout, 124) shape(weighted_jdos)
        write (stdout, 124) jdos_nbins, nspins, N_geom
124     format(3(1x, I4))
        write (stdout, '(9999(es15.8))') (((weighted_jdos(jdos_bin, N_spin, N2), N2=1, N_geom), N_spin=1, nspins) &
                                          , jdos_bin=1, jdos_nbins)
        write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
      end if

      if (optics_intraband) then
        allocate (dos_matrix_weights(size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - allocation of dos_matrix_weights failed')
        allocate (dos_at_e(3, nspins), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of dos_at_e failed')
        allocate (weighted_dos_at_e(nspins, size(matrix_weights, 5)), stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics  - allocation of weighted_dos_at_e failed')
        dos_at_e = 0.0_dp
        weighted_dos_at_e = 0.0_dp
        do N = 1, size(matrix_weights, 5)
          do N2 = 1, nbands
            dos_matrix_weights(N, N2, :, :) = matrix_weights(N2, N2, :, :, N)
          end do
        end do
        call dos_utils_calculate_at_e(efermi, dos_at_e, dos_matrix_weights, weighted_dos_at_e)
      end if

      if (on_root) then
        if (index(devel_flag, 'print_qe_constituents') > 0) then
          write (stdout, '(1x,a36,f8.4,a34)') '+------------------------ E_Fermi = ', efermi, '---------------------------------+'
          write (stdout, '(1x,a78)') '+------------------------ Printing DOS Matrix Weights -----------------------+'
          write (stdout, 125) shape(dos_matrix_weights)
          write (stdout, 125) size(matrix_weights, 5), nbands, num_kpoints_on_node(my_node_id), nspins
125       format(4(1x, I4))
          write (stdout, '(9999(es15.8))') ((((dos_matrix_weights(n_eigen, n_eigen2, N, s), s=1, nspins), N=1, &
                                              num_kpoints_on_node(my_node_id)), n_eigen2=1, nbands), n_eigen=1, &
                                            size(matrix_weights, 5))
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          write (stdout, '(1x,a78)') '+--------------------------- Printing DOS @ Energy --------------------------+'
          write (stdout, '(9(es15.8))') ((dos_at_e(i, s), i=1, 3), s=1, nspins)
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
          write (stdout, '(1x,a78)') '+----------------------- Printing Weighted DOS @ Energy ---------------------+'
          write (stdout, '(9999(es15.8))') ((weighted_dos_at_e(s, n_eigen), s=1, nspins), n_eigen=1, size(matrix_weights, 5))
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
        end if

        ! Calculate epsilon_2
        call calc_epsilon_2(weighted_jdos, weighted_dos_at_e)

        ! Calculate epsilon_1
        call calc_epsilon_1

        ! Calculate other optical properties
        call calc_refract
        call calc_absorp
        call calc_reflect

        do energy = 1, number_energies
          absorp_photo(atom, energy) = absorp(index_energy(energy))
          reflect_photo(atom, energy) = reflect(index_energy(energy))
        end do

        if (index(devel_flag, 'print_qe_constituents') > 0) then
          write (stdout, '(1x,a78)') '+-------------------- Printing Material Optical Properties ------------------+'
          write (stdout, '(1x,a78)') '+--------------------------- Printing Epsilon Array -------------------------+'
          write (stdout, 125) shape(epsilon)
          if (.not. optics_intraband) then
            write (stdout, '(9999(E17.8E3))') (((epsilon(jdos_bin, N, N2, 1), jdos_bin=1, jdos_nbins), N=1, 2), N2=1, N_geom)
          else
            write (stdout, '(9999(E17.8E3))') ((((epsilon(jdos_bin, N, N2, i), jdos_bin=1, jdos_nbins), N=1, 2), N2=1, N_geom),&
                & i=1, 3)
          end if
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

          write (stdout, '(1x,a78)') '+----------------------------- Printing Absorption --------------------------+'
          write (stdout, '(99(E17.8E3))') (absorp_photo(atom, energy), energy=1, number_energies)
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'

          write (stdout, '(1x,a78)') '+----------------------------- Printing Reflection --------------------------+'
          write (stdout, '(99(E17.8E3))') (reflect_photo(atom, energy), energy=1, number_energies)
          write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
        end if
        ! Deallocate extra arrays produced in the case of using optics_intraband
        deallocate (epsilon, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate epsilon')
        deallocate (refract, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate refract')
        deallocate (absorp, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate absorp')
        deallocate (reflect, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate reflect')
        if (optics_intraband) then
          deallocate (intra, stat=ierr)
          if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate intra')
        end if
      end if
      if (optics_intraband) then
        deallocate (dos_matrix_weights, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_matrix_weights')
        deallocate (dos_at_e, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate dos_at_e')
        deallocate (weighted_dos_at_e, stat=ierr)
        if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_dos_at_e')
      end if
      call jdos_deallocate
      deallocate (weighted_jdos, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate weighted_jdos')
      ! deallocate (E, stat=ierr)
      ! if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate E')
    end do                                        ! Loop over atoms

    call comms_bcast(absorp_photo(1, 1), max_atoms*number_energies)
    call comms_bcast(reflect_photo(1, 1), max_atoms*number_energies)

    ! Deallocating this out of the loop to reduce memory operations - could lead to higher memory consumption
    deallocate (projected_matrix_weights, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate projected_matrix_weights')

    if (index(photo_model, '1step') > 0) then
      deallocate (matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate matrix_weights')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a52,7x,f11.3,a8)') '+ Time to calculate Photoemission Optical Properties', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_photo_optics

  !***************************************************************
  subroutine calc_absorp_layer
    !***************************************************************
    ! This subroutine calculates the absorption coefficient

    use od_cell, only: atoms_pos_cart_photo
    use od_jdos_utils, only: jdos_nbins
    use od_parameters, only: devel_flag
    use od_io, only: stdout, io_error
    use od_comms, only: on_root
    implicit none
    real(kind=dp) :: I_0
    integer :: atom, i, ierr, first_atom_second_l, last_atom_secondlast_l, num_layer

    allocate (thickness_atom(max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of thickness_atom failed')
    thickness_atom = 0.0_dp

    allocate (I_layer(max_layer, number_energies), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_absorp_layer - allocation of I_layer failed')
    I_layer = 0.0_dp

    ! Calculate the thickness of each of the layers
    if (max_layer .lt. 2) then
      ! TODO: Check if this is a reasonable estimate for the single layer case
      thickness_atom = 1.5_dp
    else
      ! Setting thickness_atom for the first atom in layer # 1 and finding the first atom in the second layer
      do atom = 2, max_atoms
        if (layer(atom) .gt. 1) then
          thickness_atom(1) = ((atoms_pos_cart_photo(3, atom_order(1)) - atoms_pos_cart_photo(3, atom_order(atom)))/2)*2
          first_atom_second_l = atom
          exit
        end if
      end do
      ! Setting thickness_atom for the rest of the atoms in the first layer
      do i = 2, first_atom_second_l - 1
        thickness_atom(i) = ((atoms_pos_cart_photo(3, atom_order(i)) - &
                              atoms_pos_cart_photo(3, atom_order(first_atom_second_l)))/2)*2
      end do
      ! Setting thickness_atom for the last atom in the last layer and finding the last atom in the second to last layer
      do i = 1, max_atoms
        if (layer(max_atoms - i) .lt. layer(max_atoms)) then
          thickness_atom(max_atoms) = (abs(atoms_pos_cart_photo(3, atom_order(max_atoms)) - &
                                           atoms_pos_cart_photo(3, atom_order(max_atoms - i)))/2)*2
          last_atom_secondlast_l = max_atoms - i
          exit
        end if
      end do
      ! Setting thickness_atom for the atoms in the last layer, but not for the last atom
      do i = last_atom_secondlast_l + 1, max_atoms - 1
        thickness_atom(i) = (abs(atoms_pos_cart_photo(3, atom_order(i)) - &
                                 atoms_pos_cart_photo(3, atom_order(last_atom_secondlast_l)))/2)*2
      end do
      ! Setting thickness_atom for the atoms in between the first and last layers
      ! Formula -> abs( (z(last atom in n-1th layer) - z(atom) )/2 ) + abs( (z(first atom in n+1th layer) - z(atom) )/2 )
      do atom = first_atom_second_l, last_atom_secondlast_l
        thickness_atom(atom) = abs((atoms_pos_cart_photo(3, atom_order(sum(atoms_per_layer(1:layer(atom) - 1)))) &
                                    - atoms_pos_cart_photo(3, atom_order(atom)))/2) + &
                               abs((atoms_pos_cart_photo(3, atom_order(atom)) - &
                                    atoms_pos_cart_photo(3, atom_order(sum(atoms_per_layer(1:layer(atom))) + 1)))/2)
      end do
    end if

    I_0 = 1.0_dp
    I_layer = 1.0_dp
    do i = 1, number_energies
      I_layer(1, i) = I_0 - reflect_photo(1, i)
    end do
    if (max_layer .gt. 1) then
      do atom = first_atom_second_l, max_atoms
        do i = 1, number_energies
          ! Thickness_atom is the thickness of each layer representing the path length the light takes through the layer
          I_layer(layer(atom), i) = I_layer(layer(atom) - 1, i)* &
                                    exp(-(absorp_photo(atom, i)*thickness_atom(atom)*1E-10))
          if (I_layer(layer(atom), i) .lt. 0.0_dp) then
            I_layer(layer(atom), i) = 0.0_dp
          end if
        end do
      end do
    end if

    if (allocated(reflect_photo)) then
      deallocate (reflect_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate reflect_photo')
    end if

    ! if (allocated(atoms_per_layer)) then
    !   deallocate (atoms_per_layer, stat=ierr)
    !   if (ierr /= 0) call io_error('Error: calc_absorp_layer - failed to deallocate atoms_per_layer')
    ! end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------- Printing Intensity per Layer -----------------------+'
      write (stdout, '(1x,I4,1x,I4,1x)') jdos_nbins, max_layer
      write (stdout, '(9999(es15.8))') ((I_layer(num_layer, i), num_layer=1, max_layer), i=1, number_energies)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

  end subroutine calc_absorp_layer

  !***************************************************************
  subroutine effect_wf
    !***************************************************************
    !photo_elec_field given in eV/A

    use od_parameters, only: photo_work_function, photo_elec_field
    use od_electronic, only: efermi
    use od_constants, only: pi, epsilon_zero
    implicit none
    !  real(kind=dp) :: z

    !  z=sqrt((1/(16*pi*epsilon_zero*1E-4))/photo_elec_field)

    !  work_function_eff = photo_work_function - photo_elec_field*z -(1/(16*pi*epsilon_zero*1E-4))/z

    !  evacuum_eff = work_function_eff + efermi

    work_function_eff = photo_work_function - sqrt(photo_elec_field/(4*pi*epsilon_zero*1E-4))

    evacuum_eff = work_function_eff + efermi

  end subroutine effect_wf

  !***************************************************************
  subroutine calc_field_emission
    !***************************************************************
    ! This subroutine calculates the Schottky effect
    !photo_elec_field given in V/m

    use od_cell, only: num_kpoints_on_node
    use od_parameters, only: photo_work_function, photo_elec_field, photo_temperature, photo_surface_area
    use od_electronic, only: efermi, band_energy, nbands, nspins
    use od_io, only: io_error
    use od_comms, only: my_node_id, comms_reduce
    use od_constants, only: pi, epsilon_zero, kB, e_charge, b_factor, p1, p2, p3, p4, q1, q2, q3, q4
    implicit none
    integer :: ierr
    real(kind=dp), allocatable, dimension(:, :, :) :: field_energy
    !real(kind=dp), allocatable, dimension(:, :, :) :: tunnel_prob
    real(kind=dp), allocatable, dimension(:, :, :) :: G
    real(kind=dp), allocatable, dimension(:, :, :) :: temp_emission
    real(kind=dp) :: fermi_dirac, barrier_height, argument, exponent
    real(kind=dp) :: l_prime, p_term, q_term, v_function, transmission_prob
    integer :: N, N_spin, n_eigen

    allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of field_emission failed')
    field_emission = 0.0_dp

    allocate (field_energy(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of field_energy failed')
    field_energy = 0.0_dp

    allocate (G(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of G failed')
    G = 0.0_dp

    allocate (temp_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - allocation of temp_emission failed')
    temp_emission = 0.0_dp

    evacuum = efermi + photo_work_function

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          barrier_height = photo_work_function - (band_energy(n_eigen, N_spin, N) - efermi)
          field_energy(n_eigen, N_spin, N) = abs(evacuum - band_energy(n_eigen, N_spin, N))
          argument = (band_energy(n_eigen, N_spin, N) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a exponential value of ~1E+/-250
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 575.0_dp) then
            fermi_dirac = 0.0_dp
          elseif (argument .lt. -575.0_dp) then
            fermi_dirac = 1.0_dp
          else
            fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
          end if

          if ((photo_elec_field/(4.0_dp*pi*epsilon_zero*1E-4_dp)*photo_elec_field) .lt. (field_energy(n_eigen, N_spin, N)**2)) then
            if (barrier_height .le. 0.0_dp) then
              field_emission(n_eigen, N_spin, N) = 1.0_dp
            else
              l_prime = (e_charge/4*pi*epsilon_zero*1E-4_dp)*photo_elec_field*(1/barrier_height**2)
              p_term = 1.0_dp + (p1*l_prime) + (p2*l_prime**2.0_dp) + (p3*l_prime**3.0_dp) + (p4*l_prime**4.0_dp)
              q_term = q1 + (q2*l_prime) + (q3*l_prime**2.0_dp) + (q4*l_prime**3.0_dp)
              v_function = (1.0_dp - l_prime)*p_term + q_term*l_prime*log(l_prime)

              exponent = -1.0_dp*v_function*b_factor*(barrier_height**(2.0_dp/3.0_dp))*(1.0_dp/photo_elec_field)
              if (exponent .lt. -575.0_dp) then
                transmission_prob = 0.0_dp
              else
                transmission_prob = exp(exponent)
              end if
              field_emission(n_eigen, N_spin, N) = transmission_prob
            end if
          end if
          temp_emission(n_eigen, N_spin, N) = field_emission(n_eigen, N_spin, N)*fermi_dirac
        end do
      end do
    end do

    total_field_emission = sum(temp_emission(1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id)))/photo_surface_area
    call comms_reduce(total_field_emission, 1, "SUM")

    deallocate (field_energy, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate field_energy')

    deallocate (G, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate G')

    deallocate (temp_emission, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_field_emission - failed to deallocate temp_emission')

    !ADD COMMENT
  end subroutine calc_field_emission

  !===============================================================================
  subroutine calc_angle
    !===============================================================================
    ! This subroutine calculates the photoemission angles theta and phi
    ! Theta: angle between the photoemitted electron and the perpendicular
    !        of the surface
    ! Phi: angle between the x and y components parallel to the surface
    ! Victor Chang, 7th February 2020
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart, kpoint_r_cart
    use od_electronic, only: nbands, nspins, band_energy, band_gradient, elec_read_band_gradient, elec_read_band_curvature, &
    & band_curvature
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_momentum, devel_flag, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout, io_time
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: hbar, ev_to_j, j_to_ev, e_mass, rad_to_deg
    implicit none
    integer :: N, N_spin, n_eigen, ierr

    real(kind=dp), allocatable, dimension(:, :, :):: E_x
    real(kind=dp), allocatable, dimension(:, :, :):: E_y
    real(kind=dp) :: tol = 1.0E-10_dp
    real(kind=dp) :: time0, time1

    time0 = io_time()

    if (.not. allocated(E_transverse)) then
      allocate (E_transverse(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_transverse failed')
    end if
    E_transverse = 0.0_dp

    if (.not. allocated(theta_arpes)) then
      allocate (theta_arpes(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of theta_arpes failed')
    end if
    theta_arpes = 0.0_dp

    if (.not. allocated(phi_arpes)) then
      allocate (phi_arpes(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - allocation of phi_arpes failed')
    end if
    phi_arpes = 0.0_dp

    allocate (E_kinetic(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_kinetic failed')
    E_kinetic = 0.0_dp

    allocate (E_x(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_x failed')
    E_x = 0.0_dp

    allocate (E_y(nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - allocation of E_y failed')
    E_y = 0.0_dp

    if (index(photo_momentum, 'kp') > 0) then
      call elec_read_band_gradient
      call elec_read_band_curvature
    end if
    if (index(photo_momentum, 'operator') > 0) then
      call elec_read_band_gradient
    end if

    if (index(photo_momentum, 'crystal') > 0) call cell_calc_kpoint_r_cart

    if ((index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) .or. (index(devel_flag, 'print_qe_matrix_full') > 0 .and.&
      & on_root) .or. (index(devel_flag, 'print_qe_matrix_reduced') > 0 .and. on_root)) then
      call cell_calc_kpoint_r_cart
      write (stdout, '(a78)') "+---------------- Printing K-Points in Cartesian Coordinates ----------------+"
      do N = 1, num_kpoints_on_node(my_node_id)
        write (stdout, '(3(1x,E22.15))') kpoint_r_cart(:, N)
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          if (index(photo_momentum, 'kp') > 0) then
            E_x(n_eigen, N, N_spin) = abs &
                                      (0.5_dp*(1/(band_curvature(n_eigen, 1, 1, N, N_spin)*ev_to_j*1E-20/(hbar**2)))* &
                                       (band_gradient(n_eigen, 1, N, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
            E_y(n_eigen, N, N_spin) = abs &
                                      (0.5_dp*(1/(band_curvature(n_eigen, 2, 2, N, N_spin)*ev_to_j*1E-20/(hbar**2)))* &
                                       (band_gradient(n_eigen, 2, N, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
          end if
          if (index(photo_momentum, 'crystal') > 0) then
            E_x(n_eigen, N, N_spin) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(1, N)*1E+10)**2))*j_to_ev
            E_y(n_eigen, N, N_spin) = (((hbar**2)/(2*e_mass))*((kpoint_r_cart(2, N)*1E+10)**2))*j_to_ev
          end if
          if (index(photo_momentum, 'operator') > 0) then
            E_x(n_eigen, N, N_spin) = abs &
                                      (0.5_dp*e_mass* &
                                       (band_gradient(n_eigen, 1, N, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
            E_y(n_eigen, N, N_spin) = abs &
                                      (0.5_dp*e_mass* &
                                       (band_gradient(n_eigen, 2, N, N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
          end if
          E_transverse(n_eigen, N, N_spin) = E_x(n_eigen, N, N_spin) + E_y(n_eigen, N, N_spin)
        end do
      end do
    end do

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          if ((abs(E_x(n_eigen, N, N_spin)) .lt. tol) .and. (abs(E_y(n_eigen, N, N_spin)) .lt. tol)) then
            phi_arpes(n_eigen, N, N_spin) = 0.0_dp
          elseif ((abs(E_y(n_eigen, N, N_spin)) .lt. tol)) then
            phi_arpes(n_eigen, N, N_spin) = 90.0_dp
          else
            phi_arpes(n_eigen, N, N_spin) = atan(E_x(n_eigen, N, N_spin)/E_y(n_eigen, N, N_spin))*rad_to_deg
          end if
        end do
      end do
    end do

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands

          E_kinetic(n_eigen, N, N_spin) = (band_energy(n_eigen, N_spin, N) + temp_photon_energy - evacuum_eff)

          !! The total kinetic enery E_kinetic_total is composed of E|| and E_transvers, therefore the angle argument always
          !! has to be < 1, because theta = acos(E||/E_kinetic_total). The previous formula was: (E_kinetic(n_eigen, N, N_spin) -
          !! E_transverse(n_eigen, N, N_spin))/E_kinetic(n_eigen, N, N_spin) If the total kinetic energy is negative and
          !! E_transverse is positive, this causes the acos to be undefined, as the previous formula did not include the abs
          !!statements. Also the total kinetic energy has to be > 0 to have a physical emission
          if (E_kinetic(n_eigen, N, N_spin) .gt. 0.0_dp .and. E_kinetic(n_eigen, N, N_spin) .gt. &
          & E_transverse(n_eigen, N, N_spin)) then
            theta_arpes(n_eigen, N, N_spin) = (acos((E_kinetic(n_eigen, N, N_spin) - E_transverse(n_eigen, N, N_spin))/&
            &abs(E_kinetic(n_eigen, N, N_spin))))*rad_to_deg
          else
            theta_arpes(n_eigen, N, N_spin) = acos(0.0_dp)
          end if
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+------------------------ Printing Transverse Energy ------------------------+'
      write (stdout, '(3(1x,I4))') shape(E_transverse)
      write (stdout, '(3(1x,I4))') nbands, num_kpoints_on_node(my_node_id), nspins
      write (stdout, '(9999(es15.8))') (((E_transverse(n_eigen, N, N_spin), N_spin=1, nspins), N=1, &
                                         num_kpoints_on_node(my_node_id)), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    deallocate (E_y, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_y')

    deallocate (E_x, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_x')

    deallocate (E_kinetic, stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate E_kinetic')

    if (allocated(band_curvature)) then
      deallocate (band_curvature, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate band_curvature')
    end if

    if (allocated(band_gradient)) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate band_gradient')
    end if

    if (allocated(kpoint_r_cart)) then
      deallocate (kpoint_r_cart, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_angle - failed to deallocate kpoint_r_cart')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate Photoemission Angle', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_angle

  !***************************************************************
  subroutine calc_electron_esc
    !***************************************************************
    ! This subroutine calculates the electron escape depth

    use od_constants, only: dp, deg_to_rad
    use od_electronic, only: nbands, nspins
    use od_cell, only: num_kpoints_on_node, atoms_pos_cart_photo, atoms_label_tmp
    use od_io, only: io_error, stdout, io_time
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_imfp_const, devel_flag, iprint
    implicit none
    integer :: atom, N, N_spin, n_eigen, ierr, i
    real(kind=dp) :: exponent, time0, time1

    time0 = io_time()
    allocate (new_atoms_coordinates(3, max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of new_atoms_coordinates failed')

    !Redefine new z coordinates where the first layer is at z=0
    new_atoms_coordinates = atoms_pos_cart_photo
    do atom = 1, max_atoms
      new_atoms_coordinates(3, atom_order(atom)) = atoms_pos_cart_photo(3, atom_order(atom)) - &
                                                   (atoms_pos_cart_photo(3, atom_order(1)))
    end do

    if (.not. allocated(electron_esc)) allocate (electron_esc(nbands, nspins, &
    &num_kpoints_on_node(my_node_id), max_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    electron_esc = 0.0_dp

    if (size(photo_imfp_const, 1) .gt. 1) then
      allocate (thickness_layer(max_layer), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_layers - allocation of thickness_layer failed')
      thickness_layer = 0.0_dp

      allocate (atom_imfp(max_atoms), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_electron_esc_list - allocation of new_atoms_coordinates failed')
      atom_imfp = 0.0_dp

      ! Calculate the mean thickness of the atoms in a layer
      do atom = 1, max_atoms
        thickness_layer(layer(atom)) = thickness_layer(layer(atom)) + thickness_atom(atom)
      end do
      do i = 1, max_layer
        thickness_layer(i) = thickness_layer(i)/atoms_per_layer(i)
      end do
      write (stdout, '(1x,a78)') '+--------------- User Supplied and Calculated IMFP Constants ----------------+'
      write (stdout, '(1x,a78)') '| Atom | Atom Order | Layer | Layer Thickness | User Input IMFP | Calc. IMFP |'

      ! Calculate the layer dependent imfp constant as a list
      do atom = 1, max_atoms
        do i = 1, layer(atom)
          atom_imfp(atom) = atom_imfp(atom) + thickness_layer(i)*photo_imfp_const(i)
        end do
        atom_imfp(atom) = atom_imfp(atom)/sum(thickness_layer(1:i - 1))
        write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
          layer(atom), thickness_layer(layer(atom)), photo_imfp_const(layer(atom)), atom_imfp(atom), "    |"
225     format(1x, a1, a4, 6x, I3, 8x, I3, 6x, E13.6E3, 4x, F11.4, 3x, F11.4, a5)
      end do
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'

      do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            do atom = 1, max_atoms
              if (cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad) .gt. 0.0_dp) then
                exponent = (new_atoms_coordinates(3, atom_order(atom))/ &
                &cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad))/atom_imfp(atom)
                if (exponent .gt. -575.0_dp) then
                  electron_esc(n_eigen, N_spin, N, atom) = exp(exponent)
                else
                  electron_esc(n_eigen, N_spin, N, atom) = 0.0_dp
                end if
              end if
            end do
          end do
        end do
      end do
      deallocate (thickness_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: thickness_layer - failed to deallocate calc_elec_esc')
    elseif (size(photo_imfp_const, 1) .eq. 1) then
      do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            do atom = 1, max_atoms
              if (cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad) .gt. 0.0_dp) then
                exponent = (new_atoms_coordinates(3, atom_order(atom))/ &
                &cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad))/photo_imfp_const(1)
                if (exponent .gt. -575.0_dp) then
                  electron_esc(n_eigen, N_spin, N, atom) = exp(exponent)
                else
                  electron_esc(n_eigen, N_spin, N, atom) = 0.0_dp
                end if
              end if
            end do
          end do
        end do
      end do
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------- Printing P(Escape) per Layer -----------------------+'
      write (stdout, 125) shape(electron_esc)
      write (stdout, 125) nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms
125   format(4(1x, I4))
      write (stdout, '(9999(es15.8))') ((((electron_esc(n_eigen, N_spin, N, atom), atom=1, max_atoms), N=1, &
                                          & num_kpoints_on_node(my_node_id)), N_spin=1, nspins), n_eigen=1, nbands)
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a40,19x,f11.3,a8)') '+ Time to calculate Photoemission Escape', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_electron_esc

  !***************************************************************
  subroutine bulk_emission
    !***************************************************************
    ! This subroutine calculates the contribution from the approximated bulk material

    use od_constants, only: dp, deg_to_rad
    use od_electronic, only: nbands, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_imfp_const, photo_bulk_length, iprint
    use od_io, only: io_error, io_time, stdout
    implicit none

    real(kind=dp), dimension(:), allocatable :: bulk_light_tmp
    real(kind=dp), dimension(:, :, :, :), allocatable :: bulk_prob_tmp
    integer :: N, N_spin, n_eigen, i, num_layers, ierr
    real(kind=dp) :: exponent, time0, time1

    time0 = io_time()
    if (size(photo_imfp_const, 1) .gt. 1) then
      num_layers = int((atom_imfp(max_atoms)*photo_bulk_length)/thickness_atom(max_atoms))
    elseif (size(photo_imfp_const, 1) .eq. 1) then
      num_layers = int((photo_imfp_const(1)*photo_bulk_length)/thickness_atom(max_atoms))
    end if

    allocate (bulk_light_tmp(num_layers), stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_light_tmp failed')
    bulk_light_tmp = 0.0_dp

    allocate (bulk_prob_tmp(nbands, nspins, num_kpoints_on_node(my_node_id), num_layers), stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_prob_tmp failed')
    bulk_prob_tmp = 0.0_dp

    if (.not. allocated(bulk_prob)) allocate (bulk_prob(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - allocation of bulk_prob failed')
    bulk_prob = 0.0_dp

    if (size(photo_imfp_const, 1) .gt. 1) then
      do i = 1, num_layers
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              if (cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad) .gt. 0.0_dp) then
                exponent = (new_atoms_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms)/ &
                            cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad))/atom_imfp(max_atoms)
                ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                ! As exp(-575) is ~1E-250, this should be more than enough precision.
                if (exponent .gt. -575.0_dp) then
                  bulk_prob_tmp(n_eigen, N_spin, N, i) = exp(exponent)
                end if
              end if
            end do
          end do
        end do
      end do
      deallocate (atom_imfp, stat=ierr)
      if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate bulk_prob_tmp')

    elseif (size(photo_imfp_const, 1) .eq. 1) then
      do i = 1, num_layers
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              if (cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad) .gt. 0.0_dp) then
                exponent = (new_atoms_coordinates(3, atom_order(max_atoms)) - i*thickness_atom(max_atoms)/ &
                            cos(theta_arpes(n_eigen, N, N_spin)*deg_to_rad))/photo_imfp_const(1)
                ! This makes sure, that exp(exponent) does not underflow the dp fp value.
                ! As exp(-575) is ~1E-250, this should be more than enough precision.
                if (exponent .gt. -575.0_dp) then
                  bulk_prob_tmp(n_eigen, N_spin, N, i) = exp(exponent)
                end if
              end if
            end do
          end do
        end do
      end do

    end if

    bulk_light_tmp(1) = I_layer(layer(max_atoms), current_index)* &
                        exp(-(absorp_photo(max_atoms, current_index)*thickness_atom(max_atoms)*1E-10))
    do i = 2, num_layers
      bulk_light_tmp(i) = bulk_light_tmp(i - 1)* &
                          exp(-(absorp_photo(max_atoms, current_index)*i*thickness_atom(max_atoms)*1E-10))
    end do
    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          do i = 1, num_layers
            bulk_prob_tmp(n_eigen, N_spin, N, i) = bulk_prob_tmp(n_eigen, N_spin, N, i)*bulk_light_tmp(i)
          end do
          bulk_prob(n_eigen, N_spin, N) = sum(bulk_prob_tmp(n_eigen, N_spin, N, 1:num_layers))
        end do
      end do
    end do

    deallocate (bulk_light_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate bulk_light_tmp')

    deallocate (bulk_prob_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate bulk_prob_tmp')

    deallocate (new_atoms_coordinates, stat=ierr)
    if (ierr /= 0) call io_error('Error: bulk_emission - failed to deallocate new_atoms_coordinates')

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a38,21x,f11.3,a8)') '+ Time to calculate Bulk Photoemission', time1 - time0, ' (sec) +'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine bulk_emission

  !===============================================================================
  subroutine calc_three_step_model
    !===============================================================================
    ! This subroutine calculates the QE using the thre step model.
    ! Victor Chang, 7th February 2020
    ! edited by Felix Mildner, 03/2023
    !===============================================================================

    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient, &
      elec_read_band_curvature
    use od_comms, only: my_node_id, on_root
    use od_parameters, only: photo_surface_area, scissor_op, photo_temperature, devel_flag, photo_photon_sweep, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, io_time
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB
    implicit none
    real(kind=dp), allocatable, dimension(:, :, :, :) :: delta_temp
    real(kind=dp) :: width, norm_vac, vac_g, transverse_g, fermi_dirac, qe_factor, argument, time0, time1
    integer :: N, N2, N_spin, n_eigen, n_eigen2, atom, ierr, i

    width = (1.0_dp/11604.45_dp)*photo_temperature
    qe_factor = 1.0_dp/(2*pi*photo_surface_area)
    norm_vac = gaussian(0.0_dp, width, 0.0_dp)

    time0 = io_time()

    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of field_emission failed')
      field_emission = 0.0_dp
    end if

    if (.not. allocated(qe_tsm)) then
      allocate (qe_tsm(nbands, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - allocation of qe_tsm failed')
    end if
    qe_tsm = 0.0_dp

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+----------------- Printing Matrix Weights in 3Step Function ----------------+'
      write (stdout, '(5(1x,I4))') shape(matrix_weights)
      write (stdout, '(5(1x,I4))') nbands, nbands, num_kpoints_on_node(my_node_id), nspins, N_geom
      do N2 = 1, N_geom
        do N_spin = 1, nspins
          do N = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(es15.8))') ((matrix_weights(n_eigen, n_eigen2, N, N_spin, N2), &
            &n_eigen2=1, nbands), n_eigen=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    call photo_calculate_delta(delta_temp)

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 3Step QE ---------------------------+'
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+---------------------- Printing Delta Function Values ----------------------+'
      write (stdout, '(5(1x,I4))') shape(delta_temp)
      write (stdout, '(5(1x,I4))') nbands, nbands, num_kpoints_on_node(my_node_id), nspins
      do N_spin = 1, nspins
        do N = 1, num_kpoints_on_node(my_node_id)
          write (stdout, '(99999(es15.8))') ((delta_temp(n_eigen, n_eigen2, N_spin, N), n_eigen2=1, nbands), n_eigen=1, nbands)
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      i = 16 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
      write (stdout, '(1x,a78)') '+------------ Printing list of values going into 3step QE Values ------------+'
      write (stdout, '(1x,a261)') 'calced_qe_value - initial_state_energy - final_state_energy - matrix_weights - delta_temp -&
      & electron_esc - electrons_per_state - kpoint_weight - I_layer - qe_factor - transverse_g - vac_g - fermi_dirac -&
      & pdos_weights_atoms - pdos_weights_k_band - field_emission'
      write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', max_atoms, nbands, nbands, nspins, num_kpoints_on_node(my_node_id), i
    end if

    do atom = 1, max_atoms
      if (iprint > 2 .and. on_root) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', &
            &"Calculating atom ", atom, " of", max_atoms, "<-- QE-3S |"
      end if
      do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen2 = min_index_unocc(N_spin, N), nbands
            do n_eigen = 1, nbands

              argument = (band_energy(n_eigen, N_spin, N) - efermi)/(kB*photo_temperature)
              ! This is a bit of an arbitrary condition, but it turns out
              ! that this corresponds to a an exponent value of ~1E+/-250
              ! and this cutoff condition saves us from running into arithmetic
              ! issues when computing fermi_dirac due to possible underflow.
              if (argument .gt. 575.0_dp) then
                fermi_dirac = 0.0_dp
                exit
              elseif (argument .lt. -575.0_dp) then
                fermi_dirac = 1.0_dp
              else
                fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
              end if

              if ((temp_photon_energy - E_transverse(n_eigen, N, N_spin)) .le. (evacuum_eff - efermi)) then
                transverse_g = gaussian((temp_photon_energy - E_transverse(n_eigen, N, N_spin)), &
                                        width, (evacuum_eff - efermi))/norm_vac
              else
                transverse_g = 1.0_dp
              end if
              if ((band_energy(n_eigen, N_spin, N) + temp_photon_energy) .lt. evacuum_eff) then
                vac_g = gaussian((band_energy(n_eigen, N_spin, N) + temp_photon_energy) + &
                                 scissor_op, width, evacuum_eff)/norm_vac
              else
                vac_g = 1.0_dp
              end if

              !! this could be checked if it has an impact on the final value
              ! if (band_energy(n_eigen2, N_spin, N) .lt. efermi) cycle

              qe_tsm(n_eigen, n_eigen2, N_spin, N, atom) = &
                (matrix_weights(n_eigen, n_eigen2, N, N_spin, 1)* &
                 delta_temp(n_eigen, n_eigen2, N_spin, N)* &
                 electron_esc(n_eigen, N_spin, N, atom)* &
                 electrons_per_state*kpoint_weight(N)* &
                 (I_layer(layer(atom), current_index))* &
                 qe_factor*transverse_g*vac_g*fermi_dirac* &
                 (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom))/ &
                  pdos_weights_k_band(n_eigen, N_spin, N)))* &
                (1.0_dp + field_emission(n_eigen, N_spin, N))
              ! if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) then
              !   write (stdout, '(5(1x,I4))') atom, n_eigen, n_eigen2, N_spin, N
              !   write (stdout, '(16(1x,E17.9E3))') qe_tsm(n_eigen, n_eigen2, N_spin, N, atom), band_energy(n_eigen, N_spin, N), &
              !     band_energy(n_eigen2, N_spin, N), matrix_weights(n_eigen, n_eigen2, N, N_spin, 1), &
              !     delta_temp(n_eigen, n_eigen2, N_spin, N), electron_esc(n_eigen, N_spin, N, atom), electrons_per_state, &
              !     kpoint_weight(N), I_layer(layer(atom), current_index), qe_factor, transverse_g, vac_g, fermi_dirac, &
              !     pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom)), pdos_weights_k_band(n_eigen, N_spin, N), &
              !     field_emission(n_eigen, N_spin, N)
              ! end if
            end do
          end do
        end do
      end do
    end do

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen2 = min_index_unocc(N_spin, N), nbands
          do n_eigen = 1, nbands

            argument = (band_energy(n_eigen, N_spin, N) - efermi)/(kB*photo_temperature)
            ! This is a bit of an arbitrary condition, but it turns out
            ! that this corresponds to a an exponent value of ~1E+/-250
            ! and this cutoff condition saves us from running into arithmetic
            ! issues when computing fermi_dirac due to possible underflow.
            if (argument .gt. 575.0_dp) then
              fermi_dirac = 0.0_dp
              exit
            elseif (argument .lt. -575.0_dp) then
              fermi_dirac = 1.0_dp
            else
              fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
            end if

            if ((temp_photon_energy - E_transverse(n_eigen, N, N_spin)) .le. (evacuum_eff - efermi)) then
              transverse_g = gaussian((temp_photon_energy - E_transverse(n_eigen, N, N_spin)), &
                                      width, (evacuum_eff - efermi))/norm_vac
            else
              transverse_g = 1.0_dp
            end if
            if ((band_energy(n_eigen, N_spin, N) + temp_photon_energy) .lt. evacuum_eff) then
              vac_g = gaussian((band_energy(n_eigen, N_spin, N) + temp_photon_energy) + &
                               scissor_op, width, evacuum_eff)/norm_vac
            else
              vac_g = 1.0_dp
            end if
            qe_tsm(n_eigen, n_eigen2, N_spin, N, max_atoms + 1) = &
              (matrix_weights(n_eigen, n_eigen2, N, N_spin, 1)* &
               delta_temp(n_eigen, n_eigen2, N_spin, N)* &
               bulk_prob(n_eigen, N_spin, N)* &
               electrons_per_state*kpoint_weight(N)* &
               qe_factor*transverse_g*vac_g*fermi_dirac* &
               (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(max_atoms))/ &
                pdos_weights_k_band(n_eigen, N_spin, N)))* &
              (1.0_dp + field_emission(n_eigen, N_spin, N))
          end do
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if (allocated(delta_temp)) then
      deallocate (delta_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_three_step_model - failed to deallocate delta_temp')
    end if

    if ((index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) .or. (index(devel_flag, 'print_qe_matrix_full') > 0 .and.&
    & on_root)) then
      write (stdout, '(1x,a78)') '+----------------------- Printing Full 3step QE Matrix ----------------------+'
      write (stdout, '(5(1x,I4))') shape(qe_tsm)
      write (stdout, '(5(1x,I4))') nbands, nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1
      do atom = 1, max_atoms + 1
        do N_spin = 1, nspins
          do N = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(ES16.8E3))') ((qe_tsm(n_eigen, n_eigen2, N_spin, N, atom), &
            &n_eigen2=1, nbands), n_eigen=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 3step Photoemission', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_three_step_model

  !===============================================================================
  subroutine photo_calculate_delta(delta_temp)
    !===============================================================================
    ! It is required to evaluate the delta funcion.
    ! Victor Chang, 7th February 2020
    !===============================================================================
    use od_parameters, only: linear, fixed, adaptive, quad, iprint
    use od_electronic, only: elec_read_band_gradient, band_gradient, efermi_set
    use od_comms, only: on_root
    use od_io, only: stdout, io_error, io_time
    ! use od_cell, only: cell_volume
    use od_dos_utils, only: dos_utils_set_efermi
    use od_jdos_utils, only: setup_energy_scale, jdos_deallocate

    implicit none

    real(kind=dp) :: time0, time1
    integer       :: ierr

    real(kind=dp), intent(out), allocatable, optional    :: delta_temp(:, :, :, :)  !I've added this

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+---------------------- Calculate JDOS DELTA FUNCTION -----------------------+'
    end if

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if (quad .or. linear .or. adaptive) then
      if (.not. allocated(band_gradient)) call elec_read_band_gradient
    end if
    !-------------------------------------------------------------------------------
    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
    if (.not. efermi_set) call dos_utils_set_efermi

    time0 = io_time()

    call setup_energy_scale(E)

    if (fixed) then
      call calculate_delta('f', delta_temp)
    end if
    if (adaptive) then
      call calculate_delta('a', delta_temp)

    end if
    if (linear) then
      call calculate_delta('l', delta_temp)
    end if

    if (quad) then
      call io_error("quadratic broadening not implemented")
    end if

    call jdos_deallocate

    if (allocated(E)) then
      deallocate (E, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_calculate_delta - failed to deallocate E')
    end if

    if (allocated(band_gradient) .and. current_index .eq. number_energies) then
      deallocate (band_gradient, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_calculate_delta - failed to deallocate band_gradient')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a34,25x,f11.3,1x,a7)') &
        '+ Time to calculate Delta Function', time1 - time0, '(sec) +'
    end if
    !-------------------------------------------------------------------------------

  end subroutine photo_calculate_delta

  !===============================================================================
  subroutine calculate_delta(delta_type, delta_temp)
    !===============================================================================
    ! This subroutine evaluates the delta function between the valence band
    ! and the conduction band using the method specified in the input.
    ! Victor Chang, 7 February 2020
    !===============================================================================
    use od_comms, only: my_node_id, on_root
    use od_cell, only: num_kpoints_on_node, kpoint_grid_dim, recip_lattice
    use od_parameters, only: adaptive_smearing, fixed_smearing, iprint, &
         & finite_bin_correction, scissor_op, hybrid_linear_grad_tol, hybrid_linear, exclude_bands, num_exclude_bands, &
         & jdos_max_energy
    use od_io, only: io_error, stdout
    use od_electronic, only: band_gradient, nbands, band_energy, nspins, efermi
    use od_jdos_utils, only: jdos_nbins
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    implicit none

    integer :: ik, is, ib, jb, i, ierr
    real(kind=dp) :: cuml, width, adaptive_smearing_temp
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4), sub_cell_length(1:3)
    real(kind=dp), save                   :: delta_bins

    character(len=1), intent(in)                      :: delta_type
    real(kind=dp), intent(inout), allocatable, optional :: delta_temp(:, :, :, :)
    logical :: linear, fixed, adaptive, force_adaptive

    linear = .false.
    fixed = .false.
    adaptive = .false.

    select case (delta_type)
    case ("l")
      linear = .true.
    case ("a")
      adaptive = .true.
    case ("f")
      fixed = .true.
    case default
      call io_error(" ERROR : unknown jdos_type in calculate_delta ")
    end select

    width = 0.0_dp
    delta_bins = jdos_max_energy/real(jdos_nbins - 1, dp)

    if (linear .or. adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:), dp)/2.0_dp
    if (adaptive .or. hybrid_linear) then
      do i = 1, 3
        sub_cell_length(i) = sqrt(recip_lattice(i, 1)**2 + recip_lattice(i, 2)**2 + recip_lattice(i, 3)**2)*step(i)
      end do
      adaptive_smearing_temp = adaptive_smearing*sum(sub_cell_length)/3.0_dp
    end if

    if (fixed) width = fixed_smearing

    if (.not. allocated(delta_temp)) then
      allocate (delta_temp(nbands, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calculate_delta - allocation of delta_temp failed')
    end if
    delta_temp = 0.0_dp

    do ik = 1, num_kpoints_on_node(my_node_id)
      if (iprint > 2 .and. on_root) then
        if (mod(real(ik, dp), 10.0_dp) == 0.0_dp) write (stdout, '(1x,a1,a38,i4,a3,i4,1x,a14,2x,a11)') '^', &
             &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id), 'on this node.', "<-- Delta |"
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

            ! ! The linear method has a special way to calculate the integrated dos
            ! ! we have to take account for this here.
            ! if (linear .and. .not. force_adaptive) then
            !   delta_temp(ib, jb, ik, is) = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), photo_photon_energy, cuml)
            ! else
            !   delta_temp(ib, jb, ik, is) = gaussian((band_energy(jb,is,ik)-band_energy(ib,is,ik))+scissor_op,width,&
            !   &photo_photon_energy)
            ! end if

            ! The linear method has a special way to calculate the integrated dos
            ! we have to take account for this here.
            if (linear .and. .not. force_adaptive) then
              delta_temp(ib, jb, is, ik) = doslin(EV(0), EV(1), EV(2), EV(3), EV(4), E(current_energy_index), cuml)
            else
              delta_temp(ib, jb, is, ik) = gaussian((band_energy(jb, is, ik) - band_energy(ib, is, ik)) + scissor_op, width,&
              &E(current_energy_index))
            end if

          end do unocc_states
        end do occ_states
      end do
    end do

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if

  end subroutine calculate_delta

  !===============================================================================
  subroutine make_foptical_weights
    !===============================================================================
    ! This subroutine calclualtes the optical matrix elements for the one step
    ! photoemission model.
    ! Victor Chang, 7th February 2020
    !===============================================================================

    use od_constants, only: dp
    use od_electronic, only: nbands, nspins, num_electrons, electrons_per_state, foptical_mat
    use od_cell, only: num_kpoints_on_node, cell_get_symmetry, num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_parameters, only: optics_geom, optics_qdir, legacy_file_format, devel_flag, photo_photon_sweep
    use od_io, only: io_error, stdout
    use od_comms, only: my_node_id, on_root
    implicit none
    complex(kind=dp), dimension(3) :: g
    real(kind=dp), dimension(3) :: qdir, qdir1, qdir2
    real(kind=dp), dimension(2) :: num_occ
    real(kind=dp) :: q_weight1, q_weight2, factor
    integer :: N, i, j, N_in, N_spin, N2, N3, n_eigen, n_eigen2, num_symm, ierr

    if (.not. legacy_file_format .and. index(devel_flag, 'old_filename') > 0) then
      num_symm = 0
      call cell_get_symmetry
    end if
    num_symm = num_crystal_symmetry_operations

    num_occ = 0.0_dp
    do N_spin = 1, nspins
      num_occ(N_spin) = num_electrons(N_spin)
    end do

    if (electrons_per_state == 2) then
      num_occ(1) = num_occ(1)/2.0_dp
    end if

    ! Can I also allocate this to fome(nbands+1, num_kpts, nspins, N_geom)?
    if (.not. allocated(foptical_matrix_weights)) then
      allocate (foptical_matrix_weights(nbands + 1, nbands + 1, num_kpoints_on_node(my_node_id), nspins, N_geom), stat=ierr)
      if (ierr /= 0) call io_error('Error: make_foptical_weights - allocation of foptical_matrix_weights failed')
    end if
    foptical_matrix_weights = 0.0_dp

    if (index(optics_geom, 'polar') > 0) then
      qdir = optics_qdir
      q_weight = ((qdir(1)**2) + (qdir(2)**2) + (qdir(3)**2))**0.5_dp
      if (q_weight < 0.001_dp) &
        call io_error("Error:  please check optics_qdir, norm close to zero")
    end if

    if (index(optics_geom, 'unpolar') > 0) then
      !TO CHANGE WHEN THE light_direction IS CORRECTED
      !optics_qdir(:)=t_cart(:)
      if (optics_qdir(3) .lt. 1E-06) then
        qdir1(1) = 0.0_dp
        qdir1(2) = 0.0_dp
        qdir1(3) = 1.0_dp
      else
        qdir1(1) = 1.0_dp
        qdir1(2) = 1.0_dp
        qdir1(3) = -(optics_qdir(1) + optics_qdir(2))/optics_qdir(3)
      end if
      qdir2(1) = (optics_qdir(2)*qdir1(3)) - (optics_qdir(3)*qdir1(2))
      qdir2(2) = (optics_qdir(3)*qdir1(1)) - (optics_qdir(1)*qdir1(3))
      qdir2(3) = (optics_qdir(1)*qdir1(2)) - (optics_qdir(2)*qdir1(1))
      q_weight1 = ((qdir1(1)**2) + (qdir1(2)**2) + (qdir1(3)**2))**0.5_dp
      q_weight2 = ((qdir2(1)**2) + (qdir2(2)**2) + (qdir2(3)**2))**0.5_dp
    end if

    N_in = 1  ! 0 = no inversion, 1 = inversion
    g = 0.0_dp

    do N = 1, num_kpoints_on_node(my_node_id)                   ! Loop over kpoints
      do N_spin = 1, nspins                                    ! Loop over spins
        do n_eigen = 1, nbands                                ! Loop over state 1
          factor = 1.0_dp/(temp_photon_energy**2)
          if (index(optics_geom, 'unpolar') > 0) then
            if (num_symm == 0) then
              g(1) = (((qdir1(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                       (qdir1(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                       (qdir1(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight1)
              g(2) = (((qdir2(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                       (qdir2(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                       (qdir2(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight2)
              foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) = &
                0.5_dp*factor*(real(g(1)*conjg(g(1)), dp) + real(g(2)*conjg(g(2)), dp))
            else ! begin unpolar symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  ! Calculating foptical_matrix_weights contribution for qdir1
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir1(j))
                    end do
                  end do
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight1)
                  foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*real(g(1)*conjg(g(1)), dp)*factor
                  g(1) = 0.0_dp
                  ! Calculating foptical_matrix_weights contribution for qdir2
                  do i = 1, 3 ! if I include an extra variable I can merge this and the last do loops
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))*(crystal_symmetry_operations(j, i, N2)*qdir2(j))
                    end do
                  end do
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight2)
                  foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) + &
                    (0.5_dp/Real((num_symm*(N_in + 1)), dp))*real(g(1)*conjg(g(1)), dp)*factor
                end do
              end do
            end if !end unpolar symmetric
          elseif (index(optics_geom, 'polar') > 0) then
            if (num_symm == 0) then
              g(1) = (((qdir(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                       (qdir(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                       (qdir(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight)
              foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) = factor*real(g(1)*conjg(g(1)), dp)
            else !begin polar symmetric
              do N2 = 1, num_symm
                do N3 = 1, 1 + N_in
                  do i = 1, 3
                    qdir(i) = 0.0_dp
                    do j = 1, 3
                      qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))* &
                                (crystal_symmetry_operations(j, i, N2)*optics_qdir(j))
                    end do
                  end do
                  g(1) = 0.0_dp
                  g(1) = (((qdir(1)*foptical_mat(n_eigen, nbands + 1, 1, N, N_spin)) + &
                           (qdir(2)*foptical_mat(n_eigen, nbands + 1, 2, N, N_spin)) + &
                           (qdir(3)*foptical_mat(n_eigen, nbands + 1, 3, N, N_spin)))/q_weight)
                  foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) = &
                    foptical_matrix_weights(n_eigen, nbands + 1, N, N_spin, N_geom) + &
                    (1.0_dp/Real((num_symm*(N_in + 1)), dp))*factor*real(g(1)*conjg(g(1)), dp)
                end do
              end do
            end if !end polar symmetric
          end if ! end photo_geom
        end do       ! Loop over state 1
      end do           ! Loop over spins
    end do               ! Loop over kpoints

    if (allocated(foptical_mat) .and. current_index .eq. number_energies) then
      deallocate (foptical_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error: make_foptical_weights - failed to deallocate foptical_mat')
    end if

    if (index(devel_flag, 'print_qe_constituents') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      write (stdout, '(1x,a78)') '+------------------------- Printing Free OM Weights -------------------------+'
      write (stdout, 126) shape(foptical_matrix_weights)
      write (stdout, 126) nbands + 1, nbands + 1, num_kpoints_on_node(my_node_id), nspins, N_geom
126   format(5(1x, I4))
      do N2 = 1, N_geom
        do N_spin = 1, nspins
          do N = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(99999(es15.8))') ((foptical_matrix_weights(n_eigen, n_eigen2, N, N_spin, N2), &
                                                n_eigen2=1, nbands + 1), n_eigen=1, nbands + 1)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

  end subroutine make_foptical_weights

  !===============================================================================
  subroutine calc_one_step_model
    !===============================================================================
    ! This subroutine calculates the QE using a one step model.
    ! Victor Chang, 7th February 2020
    !===============================================================================

    use od_cell, only: num_kpoints_on_node, kpoint_weight
    use od_electronic, only: nbands, nspins, band_energy, efermi, electrons_per_state, elec_read_band_gradient,&
    & elec_read_band_curvature
    use od_comms, only: my_node_id
    use od_parameters, only: photo_surface_area, scissor_op, photo_temperature, devel_flag, photo_photon_sweep, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_comms, only: on_root
    use od_io, only: stdout, io_error, io_file_unit, io_time
    use od_jdos_utils, only: jdos_utils_calculate
    use od_constants, only: pi, kB
    implicit none
    integer :: N, N_spin, n_eigen, n_eigen2, atom, ierr, i
    real(kind=dp) :: width, norm_vac, vac_g, transverse_g, fermi_dirac, qe_factor, argument, time0, time1

    width = (1.0_dp/11604.45_dp)*photo_temperature
    qe_factor = 1.0_dp/(2*pi*photo_surface_area)
    norm_vac = gaussian(0.0_dp, width, 0.0_dp)

    time0 = io_time()
    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+--------------------------- Calculating 1Step QE ---------------------------+'
    end if

    if (.not. allocated(field_emission)) then
      allocate (field_emission(nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of field_emission failed')
      field_emission = 0.0_dp
    end if

    if (.not. allocated(qe_osm)) then
      allocate (qe_osm(nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_one_step_model - allocation of qe_osm failed')
    end if
    qe_osm = 0.0_dp

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) then
      i = 13 ! Defines the number of columns printed in the loop - needed for reshaping the data array during postprocessing
      write (stdout, '(1x,a78)') '+------------ Printing list of values going into 1step QE Values ------------+'
      write (stdout, '(1x,a222)') 'calculated_QE  foptical_matrix_weight selectron_esc electrons_per_state kpoint_weight    &
      &  I_layer      qe_factor        transverse_g         vac_g         fermi_dirac  pdos_weights_atoms pdos_weights_k_band&
      & field_emission'
      write (stdout, '(1x,a11,6(1x,I4))') 'Array Shape', i, max_atoms, nbands, nspins, num_kpoints_on_node(my_node_id)
    end if
    do atom = 1, max_atoms
      if (iprint > 2 .and. on_root) then
        write (stdout, '(1x,a1,a38,i4,a3,i4,1x,16x,a11)') ',', "Calculating atom ", atom, " of", max_atoms, "<-- QE-1S |"
      end if
      do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin = 1, nspins                    ! Loop over spins
          do n_eigen = 1, nbands
            argument = (band_energy(n_eigen, N_spin, N) - efermi)/(kB*photo_temperature)
            ! This is a bit of an arbitrary condition, but it turns out
            ! that this corresponds to a exponential value of ~1E+/-250
            ! and this cutoff condition saves us from running into arithmetic
            ! issues when computing fermi_dirac due to possible underflow.
            if (argument .gt. 555.0_dp) then
              fermi_dirac = 0.0_dp
              exit
            elseif (argument .lt. -575.0_dp) then
              fermi_dirac = 1.0_dp
            else
              fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
            end if
            if ((temp_photon_energy - E_transverse(n_eigen, N, N_spin)) .le. (evacuum_eff - efermi)) then
              transverse_g = gaussian((temp_photon_energy - E_transverse(n_eigen, N, N_spin)), &
                                      width, (evacuum_eff - efermi))/norm_vac
            else
              transverse_g = 1.0_dp
            end if
            if ((band_energy(n_eigen, N_spin, N) + temp_photon_energy) .lt. evacuum_eff) then
              vac_g = gaussian((band_energy(n_eigen, N_spin, N) + temp_photon_energy) + &
                               scissor_op, width, evacuum_eff)/norm_vac
            else
              vac_g = 1.0_dp
            end if
            n_eigen2 = nbands + 1
            qe_osm(n_eigen, N_spin, N, atom) = &
              (foptical_matrix_weights(n_eigen, n_eigen2, N, N_spin, 1)* &
               (electron_esc(n_eigen, N_spin, N, atom))* &
               electrons_per_state*kpoint_weight(N)* &
               (I_layer(layer(atom), current_index))* &
               qe_factor*transverse_g*vac_g*fermi_dirac* &
               (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom))/ &
                pdos_weights_k_band(n_eigen, N_spin, N)))* &
              (1.0_dp + field_emission(n_eigen, N_spin, N))
            ! if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root .and. .not. photo_photon_sweep) then
            !   write (stdout, '(4(1x,I4))') atom, n_eigen, N_spin, N
            !   write (stdout, '(13(7x,E16.8E4))') qe_osm(n_eigen, N_spin, N, atom), &
            !     foptical_matrix_weights(n_eigen, n_eigen2, N, N_spin, 1), &
            !    electron_esc(n_eigen, N_spin, N, atom), electrons_per_state, kpoint_weight(N), I_layer(layer(atom),current_index),&
            !     qe_factor, transverse_g, vac_g, fermi_dirac, pdos_weights_atoms(n_eigen, N_spin, N, atom_order(atom)), &
            !     pdos_weights_k_band(n_eigen, N_spin, N), field_emission(n_eigen, N_spin, N)
            ! end if
          end do
        end do
      end do
    end do
    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          argument = (band_energy(n_eigen, N_spin, N) - efermi)/(kB*photo_temperature)
          ! This is a bit of an arbitrary condition, but it turns out
          ! that this corresponds to a exponential value of ~1E+/-250
          ! and this cutoff condition saves us from running into arithmetic
          ! issues when computing fermi_dirac due to possible underflow.
          if (argument .gt. 555.0_dp) then
            fermi_dirac = 0.0_dp
            exit
          elseif (argument .lt. -575.0_dp) then
            fermi_dirac = 1.0_dp
          else
            fermi_dirac = 1.0_dp/(exp(argument) + 1.0_dp)
          end if
          if ((temp_photon_energy - E_transverse(n_eigen, N, N_spin)) .le. (evacuum_eff - efermi)) then
            transverse_g = gaussian((temp_photon_energy - E_transverse(n_eigen, N, N_spin)), &
                                    width, (evacuum_eff - efermi))/norm_vac
          else
            transverse_g = 1.0_dp
          end if
          if ((band_energy(n_eigen, N_spin, N) + temp_photon_energy) .lt. evacuum_eff) then
            vac_g = gaussian((band_energy(n_eigen, N_spin, N) + temp_photon_energy) + &
                             scissor_op, width, evacuum_eff)/norm_vac
          else
            vac_g = 1.0_dp
          end if
          n_eigen2 = nbands + 1
          qe_osm(n_eigen, N_spin, N, max_atoms + 1) = &
            (foptical_matrix_weights(n_eigen, n_eigen2, N, N_spin, 1)* &
             bulk_prob(n_eigen, N_spin, N)* &
             electrons_per_state*kpoint_weight(N)* &
             qe_factor*transverse_g*vac_g*fermi_dirac* &
             (pdos_weights_atoms(n_eigen, N_spin, N, atom_order(max_atoms))/ &
              pdos_weights_k_band(n_eigen, N_spin, N)))* &!+&
            (1.0_dp + field_emission(n_eigen, N_spin, N))
        end do
      end do
    end do

    if (index(devel_flag, 'print_qe_formula_values') > 0 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    if ((index(devel_flag, 'print_qe_constituents') > 0 .and. on_root) .or. (index(devel_flag, 'print_qe_matrix_full') > 0&
    & .and. on_root)) then
      write (stdout, '(1x,a78)') '+------------------------- Printing 1step QE Matrix -------------------------+'
      write (stdout, 125) shape(qe_osm)
      write (stdout, 125) nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1
125   format(4(1x, I4))
      do atom = 1, max_atoms + 1
        do N_spin = 1, nspins
          do N = 1, num_kpoints_on_node(my_node_id)
            write (stdout, '(9999(ES16.8E3))') (qe_osm(n_eigen, N_spin, N, atom), n_eigen=1, nbands)
          end do
        end do
      end do
      write (stdout, '(1x,a78)') '+----------------------------- Finished Printing ----------------------------+'
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a39,20x,f11.3,a8)') '+ Time to calculate 1step Photoemission', time1 - time0, ' (sec) +'
    end if

  end subroutine calc_one_step_model

  !===============================================================================
  subroutine weighted_mean_te
    !===============================================================================
    ! This subroutine calculates the weighted arithmetic mean transverse energy
    ! sum(QE*mte)/(total QE)
    ! Victor Chang, 7 February 2020
    !===============================================================================
    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nbands, nspins, elec_read_band_gradient, elec_read_band_curvature!, band_energy, efermi
    use od_comms, only: my_node_id, on_root, comms_reduce, comms_bcast
    use od_parameters, only: photo_model, iprint
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: io_error, io_file_unit, io_time, stdout
    use od_jdos_utils, only: jdos_utils_calculate

    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: te_tsm_temp
    real(kind=dp), allocatable, dimension(:, :, :, :) :: te_osm_temp
    real(kind=dp)                                     :: time0, time1
    integer :: N, N_spin, n_eigen, atom, ierr

    time0 = io_time()

    if (iprint > 1 .and. on_root) then
      write (stdout, '(1x,a78)') '+----------------------------- Calculating MTE ------------------------------+'
    end if

    if (.not. allocated(layer_qe)) then
      allocate (layer_qe(max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of layer_qe failed')
    end if
    layer_qe = 0.0_dp

    if (index(photo_model, '3step') > 0) then

      allocate (te_tsm_temp(nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of te_tsm_temp failed')
      te_tsm_temp = 0.0_dp

      ! Try : move the atom do loop to the outermost, then sum up the qe_tsm contributions from
      ! all unoccupied final states and multiply that sum by E_transverse to put into the te...
      do atom = 1, max_atoms + 1
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, min_index_unocc(N_spin, N) - 1
              !do n_eigen2 = min_index_unocc(N_spin, N), nbands
              ! if (band_energy(n_eigen2, N_spin, N) .lt. efermi) cycle ! Skip occupied final states
              te_tsm_temp(n_eigen, N, N_spin, atom) = E_transverse(n_eigen, N, N_spin) &
                                                      *sum(qe_tsm(n_eigen, min_index_unocc(N_spin, N):nbands, N_spin, N, atom))
              !end do
            end do
          end do
        end do
        ! Calculate the qe contribution of each atom/layer
        layer_qe(atom) = sum(qe_tsm(1:nbands, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), atom))
      end do

      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(layer_qe(1), max_atoms + 1, 'SUM')
      ! Calculate the total QE
      total_qe = sum(layer_qe)

      mean_te = sum(te_tsm_temp(1:nbands, 1:num_kpoints_on_node(my_node_id), 1:nspins, 1:max_atoms + 1))
      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(mean_te, 1, 'SUM')

      if (total_qe .gt. 0.0_dp) then
        mean_te = mean_te/total_qe
      else
        mean_te = 0.0_dp
      end if

      deallocate (te_tsm_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - failed to deallocate te_tsm_temp')

    elseif (index(photo_model, '1step') > 0) then

      allocate (te_osm_temp(nbands, num_kpoints_on_node(my_node_id), nspins, max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - allocation of te_osm_temp failed')
      te_osm_temp = 0.0_dp
      do atom = 1, max_atoms + 1
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              !if(band_energy(n_eigen,N_spin,N).ge.efermi) cycle
              te_osm_temp(n_eigen, N, N_spin, atom) = &
                E_transverse(n_eigen, N, N_spin)*qe_osm(n_eigen, N_spin, N, atom)
            end do
          end do
        end do
        ! Calculate the qe contribution of each atom/layer
        layer_qe(atom) = sum(qe_osm(1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), atom))
      end do

      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(layer_qe(1), max_atoms + 1, 'SUM')
      ! Calculate the total QE
      total_qe = sum(layer_qe)
      call comms_bcast(total_qe, 1)

      ! Calculate the sum of transverse E from all the bands and k-points on node
      mean_te = sum(te_osm_temp(1:nbands, 1:num_kpoints_on_node(my_node_id), 1:nspins, 1:max_atoms + 1))
      ! Sum the data from other nodes that have more k-points stored
      call comms_reduce(mean_te, 1, 'SUM')

      if (total_qe .gt. 0.0_dp) then
        mean_te = mean_te/total_qe
      else
        mean_te = 0.0_dp
      end if

      deallocate (te_osm_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: weighted_mean_te - failed to deallocate te_osm_temp')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a23,36x,f11.3,a8)') '+ Time to calculate MTE', time1 - time0, ' (sec) +'
    end if

  end subroutine weighted_mean_te

  subroutine write_qe_data
    ! This subroutine writes the calculated Photoemission data to the output file.
    ! The contents of this routine used to be part of the weighted_mean_te, but were moved
    ! here to make the subroutine names more representative of their names.

    use od_cell, only: cell_calc_kpoint_r_cart, atoms_label_tmp
    use od_comms, only: on_root
    use od_parameters, only: photo_work_function, photo_elec_field, photo_model
    use od_dos_utils, only: doslin, doslin_sub_cell_corners
    use od_algorithms, only: gaussian
    use od_io, only: stdout, io_error, io_file_unit, stdout
    use od_jdos_utils, only: jdos_utils_calculate

    integer :: atom
    if (on_root) then
      write (stdout, '(1x,a78)') '+------------------------------ Photoemission -------------------------------+'
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, 223) '| Work Function     ', photo_work_function, &
        'eV         Photon Energy', temp_photon_energy, 'eV   |'
      write (stdout, 224) '| Effective Work Function', work_function_eff, &
        'eV         Electric Field', photo_elec_field, 'V/A  |'

      if (index(photo_model, '3step') > 0) then
        write (stdout, '(1x,a78)') '| Final State : Bloch State                                                  |'
      elseif (index(photo_model, '1step') > 0) then
        write (stdout, '(1x,a78)') '| Final State : Free Electron State                                          |'
      end if
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a78)') '| Atom |  Atom Order  |   Layer   |             Quantum Efficiency           |'

      do atom = 1, max_atoms
        write (stdout, 225) "|", trim(atoms_label_tmp(atom_order(atom))), atom_order(atom), &
          layer(atom), layer_qe(atom), "      |"
      end do
      write (stdout, 226) "| Bulk", layer_qe(max_atoms + 1), &
      &"      |"

      write (stdout, 227) '| Total Quantum Efficiency (electrons/photon):', total_qe, '   |'

      write (stdout, 228) '| Weighted Mean Transverse Energy (eV):', mean_te, '      |'

      if (photo_elec_field .gt. 0.0_dp) then
        write (stdout, 229) '| Total field emission (electrons/A^2):', total_field_emission, '      |'
      end if

      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
    end if
223 format(1x, a20, f15.4, 1x, a24, f11.4, a7)
224 format(1x, a25, f10.4, 1x, a25, f10.4, a7)
225 format(1x, a1, a4, 8x, I3, 10x, I3, 16x, E24.16E3, 2x, a7)
226 format(1x, a6, 38x, E25.16E3, 2x, a7)
227 format(1x, a46, E25.16E3, a7)
228 format(1x, a39, 7x, E25.16E3, a7)
229 format(1x, a39, 7x, E25.16E3, a7)
  end subroutine write_qe_data

  subroutine binding_energy_spread
    ! TODO: Make this work well with parallelisation!!
    ! Why do we take the fixed smearing and why do we have to apply a gaussian broadening to the qe
    ! matrix? Would it make sense to apply the photo_temperature value in eV?
    ! This subroutine applies a Gaussian broadenning to the binding energy
    ! Additionally, it takes the photoemission angles theta and phi as inputs
    ! Victor Chang, 7 February 2020

    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nbands, nspins, band_energy, efermi
    use od_parameters, only: photo_work_function, fixed_smearing, photo_model, photo_theta_lower, photo_theta_upper, &
    & photo_phi_lower, photo_phi_upper
    use od_algorithms, only: gaussian
    use od_comms, only: my_node_id, comms_reduce, comms_bcast
    use od_io, only: io_error, io_file_unit
    implicit none

    real(kind=dp), allocatable, dimension(:, :, :, :) :: binding_temp
    ! real(kind=dp), allocatable, dimension(:, :, :, :) :: qe_temp
    real(kind=dp) :: qe_temp

    real(kind=dp) :: qe_norm, total_weighted
    integer :: N, N_spin, n_eigen, atom, e_scale, ierr

    max_energy = int((temp_photon_energy - photo_work_function)*1000) + 100

    allocate (t_energy(max_energy), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_spread - allocation of t_energy failed')
    t_energy = 0.0_dp

    allocate (weighted_temp(max_energy, nbands, nspins, num_kpoints_on_node(my_node_id), max_atoms + 1), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_spread - allocation of weighted_temp failed')
    weighted_temp = 0.0_dp

    allocate (binding_temp(max_energy, nbands, nspins, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_spread - allocation of binding_temp failed')
    binding_temp = 0.0_dp

    do e_scale = 1, max_energy
      t_energy(e_scale) = real(e_scale - 1, dp)/1000
    end do

    do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
      do N_spin = 1, nspins                    ! Loop over spins
        do n_eigen = 1, nbands
          do e_scale = 1, max_energy
            binding_temp(e_scale, n_eigen, N_spin, N) = &
              gaussian((efermi - band_energy(n_eigen, N_spin, N)), fixed_smearing, t_energy(e_scale))
          end do
        end do
      end do
    end do

    if (index(photo_model, '3step') > 0) then
      do atom = 1, max_atoms + 1
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              if (theta_arpes(n_eigen, N, N_spin) .ge. photo_theta_lower .and. &
                  theta_arpes(n_eigen, N, N_spin) .le. photo_theta_upper) then
                if (phi_arpes(n_eigen, N, N_spin) .ge. photo_phi_lower .and. &
                    phi_arpes(n_eigen, N, N_spin) .le. photo_phi_upper) then
                  qe_temp = sum(qe_tsm(n_eigen, 1:nbands, N_spin, N, atom))
                  do e_scale = 1, max_energy
                    weighted_temp(e_scale, n_eigen, N_spin, N, atom) = &
                      binding_temp(e_scale, n_eigen, N_spin, N)*qe_temp
                  end do
                end if
              end if
            end do
          end do
        end do
      end do

      total_weighted = sum(weighted_temp(1:max_energy, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), 1:max_atoms + 1))
      call comms_reduce(total_weighted, 1, "SUM")
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_qe/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      call comms_bcast(qe_norm, 1)
      ! Why do we need to normalise this array?
      weighted_temp = weighted_temp*qe_norm

    elseif (index(photo_model, '1step') > 0) then
      do atom = 1, max_atoms + 1
        do N = 1, num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin = 1, nspins                    ! Loop over spins
            do n_eigen = 1, nbands
              if (theta_arpes(n_eigen, N, N_spin) .ge. photo_theta_lower .and. &
                  theta_arpes(n_eigen, N, N_spin) .le. photo_theta_upper) then
                if (phi_arpes(n_eigen, N, N_spin) .ge. photo_phi_lower .and. &
                    phi_arpes(n_eigen, N, N_spin) .le. photo_phi_upper) then
                  ! if(band_energy(n_eigen,N_spin,N).ge.efermi) cycle
                  do e_scale = 1, max_energy
                    weighted_temp(e_scale, n_eigen, N_spin, N, atom) = &
                      binding_temp(e_scale, n_eigen, N_spin, N)*qe_osm(n_eigen, N_spin, N, atom)
                  end do
                end if
              end if
            end do
          end do
        end do
      end do

      total_weighted = sum(weighted_temp(1:max_energy, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), 1:max_atoms + 1))
      call comms_reduce(total_weighted, 1, "SUM")
      if (total_weighted .gt. 0.0_dp) then
        qe_norm = total_qe/total_weighted
      else
        qe_norm = 1.0_dp
      end if
      call comms_bcast(qe_norm, 1)
      ! Why do we need to normalise this array?
      weighted_temp = weighted_temp*qe_norm
    end if

    deallocate (binding_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error: binding_energy_spread - failed to deallocate binding_temp')
  end subroutine binding_energy_spread

  !***************************************************************
  subroutine write_qe_output_files
    !***************************************************************
    ! This subroutine writes either the transverse energy or the binding energy
    ! after the Gaussian broadening has been applied.
    ! Victor Chang, 7 February 2020
    ! Felix Mildner, April 2023

    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nbands, nspins
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_reduce
    use od_io, only: io_error, seedname, io_file_unit, io_date, io_time, stdout
    use od_parameters, only: write_photo_output, photo_model, iprint, devel_flag
    implicit none
    integer :: atom, ierr, e_scale, binding_unit, matrix_unit
    integer :: N, N_spin, n_eigen, kpt_total, band_num

    real(kind=dp), allocatable, dimension(:, :) :: qe_atom
    real(kind=dp) :: time0, time1
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string

    time0 = io_time()
    if (index(write_photo_output, 'qe_matrix') > 0) then
      call cell_calc_kpoint_r_cart
      kpt_total = sum(num_kpoints_on_node(0:num_nodes - 1))
      if (num_nodes .gt. 1) then
        call write_distributed_qe_data(kpt_total)
      else
        matrix_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
          filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix_final.dat'
        else
          filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix.dat'
        end if
        open (unit=matrix_unit, action='write', file=filename)
        call io_date(cdate, ctime)
        write (matrix_unit, *) '## OptaDOS Photoemission: Printing QE Matrix on ', cdate, ' at ', ctime
        write (matrix_unit, *) '## Seedname: ', trim(seedname)
        write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
        write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
        if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
          write (matrix_unit, *) '## Writing the contributions of excitations into the !!FINAL!! states'
          !   write (matrix_unit, *) '## The written values are contributions of final states to the total'
          ! elseif (index(devel_flag, 'final') .eq. 0 .and. index(photo_model, '3step') > 0) then
          !   write (matrix_unit, *) '## Writing the sum over 3-step final states contributions'
          !   write (matrix_unit, *) '## The written values are contributions of initial states to the total'
        end if
        write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'

        ! Printing out the info on root_node

        if (index(photo_model, '3step') > 0) then
          if (index(devel_flag, 'single') > 0) then
            n_eigen = len_trim(devel_flag)
            read (devel_flag(n_eigen - 2:n_eigen), *) band_num
            write (matrix_unit, '(1x,a42,1x,I3)') '## Writing contributions into final band #', band_num
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                  & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') (qe_tsm(n_eigen, band_num, N_spin, N, atom), n_eigen=1, nbands)
                end do
              end do
            end do
          else if (index(devel_flag, 'final') > 0) then
            write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
            write (matrix_unit, *) '## at a certain k-point, spin, and atom'
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                  & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') (sum(qe_tsm(1:nbands, n_eigen, N_spin, N, atom)), n_eigen=1, nbands)
                end do
              end do
            end do
          else
            write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
            write (matrix_unit, *) '## at a certain k-point, spin, and atom'
            write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
                  & ')'
            do atom = 1, max_atoms + 1
              if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
              do N = 1, num_kpoints_on_node(my_node_id)
                do N_spin = 1, nspins
                  write (matrix_unit, '(9999(ES16.8E3))') (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N, atom)), n_eigen=1, nbands)
                end do
              end do
            end do
          end if
        elseif (index(photo_model, '1step') > 0) then
          write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
          write (matrix_unit, *) '## at a certain k-point, spin, and atom'
          write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms,&
               & ')'
          do atom = 1, max_atoms + 1
            if (atom .eq. max_atoms + 1) write (matrix_unit, *) '## Bulk Contribution:'
            do N = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') (qe_osm(n_eigen, N_spin, N, atom), n_eigen=1, nbands)
              end do
            end do
          end do
        end if
        close (unit=matrix_unit)
      end if
    end if

    if (index(write_photo_output, 'e_bind') > 0) then

      allocate (qe_atom(max_energy, max_atoms + 1), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - allocation of qe_atom failed')
      qe_atom = 0.0_dp

      do atom = 1, max_atoms + 1
        do e_scale = 1, max_energy !loop over transverse energy
          qe_atom(e_scale, atom) = &
            sum(weighted_temp(e_scale, 1:nbands, 1:nspins, 1:num_kpoints_on_node(my_node_id), atom))
        end do
      end do

      if (num_nodes .gt. 1) then
        call comms_reduce(qe_atom(1, 1), max_energy*(max_atoms + 1), "SUM")
      end if

      if (on_root) then
        binding_unit = io_file_unit()
        write (char_e, '(F7.3)') temp_photon_energy
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//&
        &'_binding_energy.dat'
        open (unit=binding_unit, action='write', file=filename)
        call io_date(cdate, ctime)
        write (binding_unit, *) '## OptaDOS Photoemission: Printing Binding Energy Spread on ', cdate, ' at ', ctime
        write (binding_unit, *) '## Seedname: ', trim(seedname)
        write (binding_unit, *) '## Photoemission Model: ', trim(photo_model)
        write (binding_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
        write (binding_unit, *) '## Transverse Energy (TE) [eV] | Total QE from sum(atoms + bulk) @ TE | Contributions from: atom1 &
        &| atom2 | ... | bulk | '

        do e_scale = 1, max_energy
          write (binding_unit, '(1x,ES13.6E2,2x,ES19.12E3,1x,999(1x,ES19.12E3))') t_energy(e_scale), &
          &sum(qe_atom(e_scale, 1:max_atoms + 1)), qe_atom(e_scale, 1:max_atoms + 1)
        end do

        close (unit=binding_unit)
      end if
    end if

    if (allocated(weighted_temp)) then
      deallocate (weighted_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate weighted_temp')
    end if

    if (allocated(qe_atom)) then
      deallocate (qe_atom, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate qe_atom')
    end if

    if (allocated(t_energy)) then
      deallocate (t_energy, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_qe_output_files - failed to deallocate t_energy')
    end if

    time1 = io_time()
    if (on_root .and. iprint > 1) then
      write (stdout, '(1x,a78)') '+----------------------------------------------------------------------------+'
      write (stdout, '(1x,a32,27x,f11.3,a8)') '+ Time to write the output files', time1 - time0, ' (sec) +'
    end if

  end subroutine write_qe_output_files

  subroutine write_distributed_qe_data(kpt_total)
    !***************************************************************
    ! This subroutine writes the distributed qe tensor to a single file.
    ! To save on required memory the output file is accessed by each MPI process in turn
    ! and writes its values/contents one after the other.
    ! F. Mildner, June 2023

    use od_cell, only: num_kpoints_on_node, cell_calc_kpoint_r_cart
    use od_electronic, only: nspins, nbands
    use od_comms, only: my_node_id, on_root, num_nodes, comms_send, comms_recv, root_id, comms_bcast
    use od_io, only: io_error, io_file_unit, io_date, io_time, seedname
    use od_parameters, only: photo_model, devel_flag

    implicit none
    real(kind=dp), dimension(:, :, :), allocatable :: qe_mat_temp
    real(kind=dp), dimension(:, :, :, :), allocatable :: tsm_reduced
    integer, intent(in)                         :: kpt_total
    character(len=99)                           :: filename
    character(len=10)                           :: char_e
    character(len=9)                            :: ctime             ! Temp. time string
    character(len=11)                           :: cdate             ! Temp. date string
    integer:: N, N_spin, n_eigen, atom, token, matrix_unit, ierr, inode

    ! On root open file and write header

    if (on_root) then
      ! Writing header to output file
      write (char_e, '(F7.3)') temp_photon_energy
      if (index(devel_flag, 'final') > 0 .and. index(photo_model, '3step') > 0) then
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix_final.dat'
      else
        filename = trim(seedname)//'_'//trim(photo_model)//'_'//trim(adjustl(char_e))//'_qe_matrix.dat'
      end if
      matrix_unit = io_file_unit()
      open (unit=matrix_unit, action='write', file=filename)
      call io_date(cdate, ctime)
      write (matrix_unit, *) '## OptaDOS Photoemission: Printing QE Matrix on ', cdate, ' at ', ctime
      write (matrix_unit, *) '## Seedname: ', trim(seedname)
      write (matrix_unit, *) '## Photoemission Model: ', trim(photo_model)
      write (matrix_unit, *) '## Photon Energy: ', trim(adjustl(char_e))
      write (matrix_unit, *) '## Find band energies and fractional k-point coordinates in: ', trim(seedname), '.bands'
      write (matrix_unit, *) '## (Reduced) QE Matrix where each row contains the contributions from each band'
      write (matrix_unit, *) '## at a certain k-point, spin, and atom'
      write (matrix_unit, '(1x,a31,4(1x,I5),1x,1a)') '## (Reduced) QE Matrix Shape: (', nbands, nspins, kpt_total, max_atoms, ')'
      allocate (qe_mat_temp(nbands, nspins, num_kpoints_on_node(0)), stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate qe_mat_temp on root')
      token = -1
    end if

    ! allocate and sum the 3step qe matrix on non-root
    if (.not. on_root) then
      if (index(photo_model, '3step') > 0) then
        allocate (tsm_reduced(nbands, nspins, num_kpoints_on_node(0), max_atoms + 1), stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to allocate tsm_reduced')
        if (index(devel_flag, 'final') > 0) then
          tsm_reduced = sum(qe_tsm, dim=1)
        else
          tsm_reduced = sum(qe_tsm, dim=2)
        end if
      end if
    end if
    ! For each atom until max_atoms+1
    do atom = 1, max_atoms + 1
      ! On non root nodes
      if (.not. on_root) then
        ! - wait for the token
        call comms_recv(token, 1, 0)
        ! - send the respective qe_matrix for that specific atom
        if (index(photo_model, '3step') > 0) then
          call comms_send(tsm_reduced(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
        elseif (index(photo_model, '1step') > 0) then
          call comms_send(qe_osm(1, 1, 1, atom), nbands*nspins*num_kpoints_on_node(my_node_id), 0)
        end if
        ! - send token back to root node
        call comms_send(token, 1, 0)
        ! On root node
      elseif (on_root) then
        do inode = 1, num_nodes - 1
          ! - send to the token to notes in turn
          call comms_send(token, 1, inode)
          ! - receive the qe_matrix from the other notes and write it to the file
          call comms_recv(qe_mat_temp(1, 1, 1), nbands*nspins*num_kpoints_on_node(inode), inode)
          ! write out the qe_matrix to the file
          do N = 1, num_kpoints_on_node(inode)
            do N_spin = 1, nspins
              write (matrix_unit, '(9999(ES16.8E3))') (qe_mat_temp(n_eigen, N_spin, N), n_eigen=1, nbands)
            end do
          end do
          ! - receive the token from a node
          call comms_recv(token, 1, inode)
        end do
        ! - write root qe_matrix elements
        if (index(photo_model, '3step') > 0) then
          if (index(devel_flag, 'final') > 0) then
            do N = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') (sum(qe_tsm(1:nbands, n_eigen, N_spin, N, atom)), n_eigen=1, nbands)
              end do
            end do
          else
            do N = 1, num_kpoints_on_node(my_node_id)
              do N_spin = 1, nspins
                write (matrix_unit, '(9999(ES16.8E3))') (sum(qe_tsm(n_eigen, 1:nbands, N_spin, N, atom)), n_eigen=1, nbands)
              end do
            end do
          end if
        elseif (index(photo_model, '1step') > 0) then
          do N = 1, num_kpoints_on_node(my_node_id)
            do N_spin = 1, nspins
              write (matrix_unit, '(9999(ES16.8E3))') (qe_osm(n_eigen, N_spin, N, atom), n_eigen=1, nbands)
            end do
          end do
        end if
        ! Write header for bulk contrib using root node
        if (atom .eq. max_atoms) write (matrix_unit, '(1x,a21)') '## Bulk Contribution:'
      end if
    end do
    if (on_root) then
      close (unit=matrix_unit)
      deallocate (qe_mat_temp, stat=ierr)
      if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to deallocate qe_mat_temp')
    elseif (.not. on_root) then
      if (index(photo_model, '3step') > 0) then
        deallocate (tsm_reduced, stat=ierr)
        if (ierr /= 0) call io_error('Error: write_distributed_qe_data - failed to deallocate tsm_reduced')
      end if
    end if
  end subroutine write_distributed_qe_data

  subroutine photo_deallocate
    !***************************************************************
    ! This subroutine deallocates all the quantities which have not
    ! been deallocated yet

    use od_io, only: io_error
    use od_electronic, only: foptical_mat
    implicit none
    integer :: ierr

    if (allocated(phi_arpes)) then
      deallocate (phi_arpes, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate phi_arpes')
    end if

    if (allocated(theta_arpes)) then
      deallocate (theta_arpes, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate theta_arpes')
    end if

    if (allocated(refract)) then
      deallocate (refract, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate refract')
    end if

    if (allocated(absorp)) then
      deallocate (absorp, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate absorp')
    end if

    if (allocated(electron_esc)) then
      deallocate (electron_esc, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate electron_esc')
    end if

    if (allocated(layer)) then
      deallocate (layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate layer')
    end if

    if (allocated(layer_qe)) then
      deallocate (layer_qe, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate layer_qe')
    end if

    if (allocated(imfp_val)) then
      deallocate (imfp_val, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate imfp_val')
    end if

    if (allocated(reflect)) then
      deallocate (reflect, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate reflect')
    end if

    if (allocated(matrix_weights)) then
      deallocate (matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: calc_photo_optics - failed to deallocate matrix_weights')
    end if

    if (allocated(E_transverse)) then
      deallocate (E_transverse, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate E_transverse')
    end if

    if (allocated(absorp_photo)) then
      deallocate (absorp_photo, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate absorp_photo')
    end if

    if (allocated(atom_order)) then
      deallocate (atom_order, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate atom_order')
    end if

    if (allocated(atoms_per_layer)) then
      deallocate (atoms_per_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate atoms_per_layer')
    end if

    if (allocated(thickness_atom)) then
      deallocate (thickness_atom, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate thickness_atom')
    end if

    if (allocated(pdos_weights_atoms)) then
      deallocate (pdos_weights_atoms, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_atoms')
    end if

    if (allocated(pdos_weights_k_band)) then
      deallocate (pdos_weights_k_band, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_k_band')
    end if

    if (allocated(index_energy)) then
      deallocate (index_energy, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate index_energy')
    end if

    if (allocated(I_layer)) then
      deallocate (I_layer, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate I_layer')
    end if

    if (allocated(field_emission)) then
      deallocate (field_emission, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate field_emission')
    end if

    if (allocated(bulk_prob)) then
      deallocate (bulk_prob, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate bulk_prob')
    end if

    if (allocated(qe_tsm)) then
      deallocate (qe_tsm, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate qe_tsm')
    end if

    if (allocated(qe_osm)) then
      deallocate (qe_osm, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate qe_osm')
    end if

    if (allocated(foptical_mat)) then
      deallocate (foptical_mat, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate foptical_mat')
    end if

    if (allocated(foptical_matrix_weights)) then
      deallocate (foptical_matrix_weights, stat=ierr)
      if (ierr /= 0) call io_error('Error: photo_deallocate - failed to deallocate foptical_matrix_weights')
    end if

  end subroutine photo_deallocate

end module od_photo
