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

module od_core

  use od_constants, only: dp
  implicit none
  private
  public :: core_calculate

  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: matrix_weights
  real(kind=dp), allocatable, public, dimension(:, :, :) :: weighted_dos
  real(kind=dp), allocatable, public, dimension(:, :, :) :: weighted_dos_broadened

contains

  subroutine core_calculate
    use od_electronic, only: elec_read_elnes_mat, efermi_set
    use od_dos_utils, only: dos_utils_calculate, dos_utils_set_efermi, &
    & dos_utils_compute_bandgap
    use od_comms, only: on_root
    use od_io, only: stdout
    use od_parameters, only: core_LAI_broadening, LAI_gaussian, LAI_lorentzian, &
    & set_efermi_zero, LAI_lorentzian_scale, core_chemical_shift

    implicit none

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                            Core Loss Calculation                           +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '|                                                                            |'
    end if

    ! read in the core matrix elements from disk
    call elec_read_elnes_mat
    !    (elnes_mat(orb,nb,indx,nk,ns),indx=1,3)

    if (.not. efermi_set) call dos_utils_set_efermi

    call core_prepare_matrix_elements

    call dos_utils_calculate(matrix_weights, weighted_dos)

    ! Lifetime and instrumental broadening
    if (core_LAI_broadening .eqv. .true.) then
      allocate (weighted_dos_broadened(size(weighted_dos, 1), size(weighted_dos, 2), size(weighted_dos, 3)))
      weighted_dos_broadened = 0.0_dp
      if (LAI_lorentzian .or. (LAI_lorentzian_scale .gt. 0.00001_dp)) call core_lorentzian
      if (LAI_gaussian) call core_gaussian
    end if

    if (set_efermi_zero .and. .not. efermi_set) call dos_utils_set_efermi
    if (on_root) then
      call write_core
    end if

  end subroutine core_calculate

  ! Private routines

  subroutine core_prepare_matrix_elements
    use od_electronic, only: elnes_mat, elnes_mwab, nbands, nspins, num_electrons, electrons_per_state, &
      efermi, band_energy
    use od_comms, only: my_node_id
    use od_cell, only: num_kpoints_on_node, cell_get_symmetry, &
      num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_parameters, only: core_geom, core_qdir, core_type, legacy_file_format, devel_flag
    use od_io, only: io_error

    real(kind=dp), dimension(3) :: qdir
    real(kind=dp) :: q_weight
    integer :: N, N_spin, n_eigen, orb, ierr, num_sym, j, N2, N3, i, N_in
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp) :: g

    num_occ = 0.0_dp

    if (.not. legacy_file_format .and. index(devel_flag, 'old_filename') > 0) then
      call cell_get_symmetry
    end if
    num_sym = num_crystal_symmetry_operations

    if (num_sym > 1) call io_error('Error: Core loss not currently able to deal with symmetry. Please &
          &re-run CASTEP without symmetry.')

    allocate (matrix_weights(elnes_mwab%norbitals, elnes_mwab%nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_prepare_matrix_elements - allocation failed for matrix_weights')
    matrix_weights = 0.0_dp

    N_in = 1  ! 0 = no inversion, 1 = inversion

    if (index(core_geom, 'polar') > 0) then
      qdir = core_qdir
      q_weight = ((qdir(1)**2.0_dp) + (qdir(2)**2.0_dp) + (qdir(3)**2.0_dp))**0.5_dp
      if (q_weight < 0.001_dp) &
        call io_error("Error: core_prepare_matrix_elements.  please check core_qdir, norm close to zero")
    end if

    do N = 1, num_kpoints_on_node(my_node_id)                      ! Loop over kpoints
      do N_spin = 1, nspins                                    ! Loop over spins
        do n_eigen = 1, nbands                                ! Loop over state 1
          if (band_energy(n_eigen, N_spin, N) < efermi .and. core_type == 'absorption') cycle ! XES
          if (band_energy(n_eigen, N_spin, N) > efermi .and. core_type == 'emission') cycle   ! ELNES / XANES
          do orb = 1, elnes_mwab%norbitals
            if (index(core_geom, 'polar') > 0) then
              if (num_sym == 0) then
                g = (((qdir(1)*elnes_mat(orb, n_eigen, 1, N, N_spin)) + &
                      (qdir(2)*elnes_mat(orb, n_eigen, 2, N, N_spin)) + &
                      (qdir(3)*elnes_mat(orb, n_eigen, 3, N, N_spin)))/q_weight)
                matrix_weights(orb, n_eigen, N, N_spin) = real(g*conjg(g), dp)
              else
                do N2 = 1, num_sym
                  do N3 = 1, 1 + N_in
                    do i = 1, 3
                      qdir(i) = 0.0_dp
                      do j = 1, 3
                        qdir(i) = qdir(i) + ((-1.0_dp)**(N3 + 1))* &
                                  (crystal_symmetry_operations(j, i, N2)*core_qdir(j))
                      end do
                    end do
                    g = 0.0_dp
                    g = (((qdir(1)*elnes_mat(orb, n_eigen, 1, N, N_spin)) + &
                          (qdir(2)*elnes_mat(orb, n_eigen, 2, N, N_spin)) + &
                          (qdir(3)*elnes_mat(orb, n_eigen, 3, N, N_spin)))/q_weight)
                    matrix_weights(orb, n_eigen, N, N_spin) = &
                      matrix_weights(orb, n_eigen, N, N_spin) + &
                      (1.0_dp/Real((num_sym*(N_in + 1)), dp))*real(g*conjg(g), dp)
                  end do
                end do
              end if
            else ! isotropic average
              matrix_weights(orb, n_eigen, N, N_spin) = real( &
                & elnes_mat(orb, n_eigen, 1, N, N_spin)*conjg(elnes_mat(orb, n_eigen, 1, N, N_spin)) + &
                & elnes_mat(orb, n_eigen, 2, N, N_spin)*conjg(elnes_mat(orb, n_eigen, 2, N, N_spin)) + &
                & elnes_mat(orb, n_eigen, 3, N, N_spin)*conjg(elnes_mat(orb, n_eigen, 3, N, N_spin)), dp) &
                & /3.0_dp
              !                matrix_weights(orb,n_eigen,N,N_spin) = real(g*conjg(g),dp)
              !           matrix_weights(n_eigen,n_eigen2,N,N_spin) = 1.0_dp  !
            end if
          end do
        end do
      end do
    end do

  end subroutine core_prepare_matrix_elements

  subroutine write_core
    !*************************************************************************
    ! This subroutine writes out the Core loss function
    !-------------------------------------------------------------------------
    ! Adapted by A F Harper to include an E_shift to account for core hole
    !=========================================================================

    use od_constants, only: bohr2ang, periodic_table_name, pi
    use od_parameters, only: dos_nbins, core_LAI_broadening, LAI_gaussian, LAI_gaussian_width, &
      LAI_lorentzian, LAI_lorentzian_scale, LAI_lorentzian_width, LAI_lorentzian_offset, output_format, &
      set_efermi_zero, core_chemical_shift
    use od_electronic, only: elnes_mwab, elnes_orbital, efermi, efermi_set, nspins
    use od_io, only: seedname, io_file_unit, io_error
    use od_dos_utils, only: E, dos_utils_set_efermi, vbm_energy, cbm_energy
    use od_cell, only: num_species, atoms_symbol, atoms_label, cell_volume
    use xmgrace_utils

    integer :: N
    real(kind=dp) ::dE, min_x, min_y, max_x, max_y, range
    integer :: core_unit, orb, ierr, loop, loop2, counter, num_edge, num_sites
    character(len=20) :: temp
    character(len=40) :: temp2
    character(len=10), allocatable :: elnes_symbol(:)
    character(len=10), allocatable :: elnes_label(:)
    integer, allocatable :: edge_shell(:), edge_am(:), edge_num_am(:), edge_list(:, :)
    integer, allocatable :: edge_species(:), edge_rank_in_species(:)
    integer, allocatable :: ion_species(:), ion_num_in_species(:)
    character(len=40), allocatable :: edge_name(:)
    logical :: found
    real(kind=dp), allocatable :: dos_temp(:, :), dos_temp2(:, :)
    real(kind=dp), allocatable :: E_shift(:)
    real(kind=dp) :: epsilon2_const
    real(kind=dp), parameter :: epsilon_0 = 8.8541878176E-12_dp
    real(kind=dp), parameter :: e_charge = 1.602176487E-19_dp

    allocate (E_shift(dos_nbins), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating E_shift in core_write')
    if (set_efermi_zero) then
      E_shift = E - efermi
    else
      E_shift = E
    end if

    if (nspins == 1) then
      allocate (dos_temp(dos_nbins, 1), stat=ierr)
    else
      allocate (dos_temp(dos_nbins, 3), stat=ierr)
    end if
    if (ierr /= 0) call io_error('Error: core_write - allocation of dos_temp failed')

    if (nspins == 1) then
      allocate (dos_temp2(dos_nbins, 1), stat=ierr)
    else
      allocate (dos_temp2(dos_nbins, 3), stat=ierr)
    end if
    if (ierr /= 0) call io_error('Error: core_write - allocation of dos_temp2 failed')

    allocate (elnes_symbol(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of elnes_symbol failed')
    allocate (elnes_label(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of elnes_label failed')

    dE = E(2) - E(1)

    ! This next bit of convoluted code attempts to figure out what order the
    ! orbitals were in the elnes file. It will be the order castep labels the atoms
    ! which is not the same as the cell file. Not that any species that are defined
    ! using labels such as B:ext B:1 etc will appear after the other elements.

    counter = 1
    do loop = 1, num_species
      do loop2 = 1, 109
        if (atoms_label(loop) == periodic_table_name(loop2)) then
          elnes_symbol(counter) = periodic_table_name(loop2)
          elnes_label(counter) = ''
          counter = counter + 1
          exit
          !check atom count here
        end if
      end do
    end do
    if (counter < num_species + 1) then
      do loop = 1, num_species
        found = .false.
        do loop2 = 1, 109
          if (atoms_label(loop) == periodic_table_name(loop2)) then
            found = .true.
            exit
          end if
        end do
        if (.not. found) then
          do loop2 = 1, 109
            if (atoms_symbol(loop) == periodic_table_name(loop2)) then
              elnes_symbol(counter) = periodic_table_name(loop2)
              elnes_label(counter) = atoms_label(loop)
              counter = counter + 1
              exit
            end if
          end do
        end if
      end do
    end if

!!$
!!$    ! Open the output file
!!$    core_unit = io_file_unit()
!!$    open(unit=core_unit,action='write',file=trim(seedname)//'_core.dat')
!!$
!!$    ! Write into the output file
!!$    write(core_unit,*)'#*********************************************'
!!$    write(core_unit,*)'#            Core loss function               '
!!$    write(core_unit,*)'#*********************************************'
!!$    write(core_unit,*)'#'
!!$    if(core_LAI_broadening) then
!!$       if(LAI_gaussian) write(core_unit,*)'# Gaussian broadening: FWHM', LAI_gaussian_width
!!$       if(LAI_lorentzian) then
!!$          write(core_unit,*)'# Lorentzian broadening included'
!!$          write(core_unit,*)'# Lorentzian scale ', LAI_lorentzian_scale
!!$          write(core_unit,*)'# Lorentzian offset ', LAI_lorentzian_offset
!!$          write(core_unit,*)'# Lorentzian width ', LAI_lorentzian_width
!!$       end if
!!$    end if
!!$    write(core_unit,*)'#'
!!$
!!$
!!$    do orb=1,elnes_mwab%norbitals   ! writing out doesn't include spin at the moment.
!!$
!!$       write(temp,'(a2,i1,a17)') 'n=',elnes_orbital%shell(orb),' ang= '//trim(elnes_orbital%am_channel_name(orb))
!!$
!!$       write(temp2,'(a5,a2,1x,i0,a28)') 'Ion: ',trim(elnes_symbol(elnes_orbital%species_no(orb))),&
!!$        & elnes_orbital%rank_in_species(orb),' State: '//trim(temp)
!!$       write(core_unit,*) '# ',trim(temp2)
!!$
!!$
!!$       do N=1,dos_nbins
!!$          if (core_LAI_broadening) then
!!$             write(core_unit,*)E(N),weighted_dos(N,1,orb),weighted_dos_broadened(N,1,orb)
!!$          else
!!$             write(core_unit,*)E(N),weighted_dos(N,1,orb)
!!$          end if
!!$       end do
!!$       write(core_unit,*)''
!!$    end do
!!$
!!$    close(core_unit)

    ! Converts units, note I don't have to worry about this in optics.f90 at electronic does it
    weighted_dos = weighted_dos*bohr2ang**2
    if (core_LAI_broadening) weighted_dos_broadened = weighted_dos_broadened*bohr2ang**2

    ! Units are currently (ang^2)(eV^-1) so need to multiply by a factor so we have a dimensionless epsilon_2
    epsilon2_const = (e_charge*pi*1E-20)/(cell_volume*1E-30*epsilon_0)
    weighted_dos = weighted_dos*epsilon2_const
    if (core_LAI_broadening) weighted_dos_broadened = weighted_dos_broadened*epsilon2_const

    allocate (ion_species(elnes_mwab%norbitals))
    allocate (ion_num_in_species(elnes_mwab%norbitals))
    ion_species = 0; ion_num_in_species = 0
    !Find out what ions we have
    do loop = 1, elnes_mwab%norbitals
      if (loop == 1) then
        ion_species(1) = elnes_orbital%species_no(1)
        ion_num_in_species(1) = elnes_orbital%rank_in_species(1)
        counter = 1
      else
        found = .false.
        do loop2 = 1, counter
          if (elnes_orbital%species_no(loop) == ion_species(loop2) .and.&
               & elnes_orbital%rank_in_species(loop) == ion_num_in_species(loop2)) then
            found = .true.
          end if
        end do
        if (.not. found) then
          counter = counter + 1
          ion_species(counter) = elnes_orbital%species_no(loop)
          ion_num_in_species(counter) = elnes_orbital%rank_in_species(loop)
        end if
      end if
    end do
    num_sites = counter

    ! We allocate these arrays as the max possible size, and just fill in the bits we need
    allocate (edge_species(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_species failed')
    allocate (edge_rank_in_species(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_rank_in_species failed')
    allocate (edge_shell(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_shell failed')
    allocate (edge_am(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_am failed')
    allocate (edge_num_am(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_num_am failed')
    allocate (edge_list(elnes_mwab%norbitals, 7), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_list failed')
    edge_species = 0; edge_rank_in_species = 0; edge_shell = 0; edge_am = 0; edge_num_am = 0; edge_list = 0

    counter = 1
    ! Find out how many edges
    do loop = 1, elnes_mwab%norbitals
      if (loop == 1) then
        edge_species(counter) = elnes_orbital%species_no(loop)
        edge_rank_in_species(counter) = elnes_orbital%rank_in_species(loop)
        edge_shell(counter) = elnes_orbital%shell(loop)
        edge_am(counter) = elnes_orbital%am_channel(loop)
        edge_num_am(counter) = edge_num_am(counter) + 1
        edge_list(counter, edge_num_am(counter)) = loop
      else
        ! else check if we have this am state
        found = .false.
        do loop2 = 1, counter
          if (edge_species(loop2) == elnes_orbital%species_no(loop) .and.&
               & edge_rank_in_species(loop2) == elnes_orbital%rank_in_species(loop) .and. &
               edge_shell(loop2) == elnes_orbital%shell(loop) .and. edge_am(loop2) == elnes_orbital%am_channel(loop)) then
            edge_num_am(counter) = edge_num_am(counter) + 1
            edge_list(counter, edge_num_am(counter)) = loop
            found = .true.
          end if
        end do
        if (.not. found) then
          counter = counter + 1
          edge_species(counter) = elnes_orbital%species_no(loop)
          edge_rank_in_species(counter) = elnes_orbital%rank_in_species(loop)
          edge_shell(counter) = elnes_orbital%shell(loop)
          edge_am(counter) = elnes_orbital%am_channel(loop)
          edge_num_am(counter) = edge_num_am(counter) + 1
          edge_list(counter, edge_num_am(counter)) = loop
        end if
      end if
    end do
    num_edge = counter
    !
    allocate (edge_name(num_edge), stat=ierr)
    if (ierr /= 0) call io_error('Error: core_write - allocation of edge_name failed')
    ! fill in edge name
    do loop = 1, num_edge
      if (edge_shell(loop) == 1) then
        temp = 'K1'
      elseif (edge_shell(loop) == 2) then
        if (edge_am(loop) == 0) then
          temp = 'L1'
        elseif (edge_am(loop) == 1) then
          temp = 'L2,3'
        end if
      elseif (edge_shell(loop) == 3) then
        if (edge_am(loop) == 0) then
          temp = 'M1'
        elseif (edge_am(loop) == 1) then
          temp = 'M2,3'
        elseif (edge_am(loop) == 2) then
          temp = 'M4,5'
        end if
      elseif (edge_shell(loop) == 4) then
        if (edge_am(loop) == 0) then
          temp = 'N1'
        elseif (edge_am(loop) == 1) then
          temp = 'N2,3'
        elseif (edge_am(loop) == 2) then
          temp = 'N4,5'
        elseif (edge_am(loop) == 3) then
          temp = 'N6,7'
        end if
      elseif (edge_shell(loop) == 5) then
        if (edge_am(loop) == 0) then
          temp = 'O1'
        elseif (edge_am(loop) == 1) then
          temp = 'O2,3'
        elseif (edge_am(loop) == 2) then  ! after this point I think we've drifted beyond what is physical!
          temp = 'O4,5'
        elseif (edge_am(loop) == 3) then
          temp = 'O6,7'
        end if
      elseif (edge_shell(loop) == 6) then
        if (edge_am(loop) == 0) then
          temp = 'P1'
        elseif (edge_am(loop) == 1) then
          temp = 'P2,3'
        elseif (edge_am(loop) == 2) then
          temp = 'P4,5'
        elseif (edge_am(loop) == 3) then
          temp = 'P6,7'
        end if
      end if

      if (elnes_label(edge_species(loop)) == '') then
        write (edge_name(loop), '(a2,1x,i0,1x,a5)') trim(elnes_symbol(edge_species(loop))), &
        & edge_rank_in_species(loop), trim(temp)
      else
        write (edge_name(loop), '(a2,1x,i0,1x,a5,a10)') trim(elnes_symbol(edge_species(loop))), &
          edge_rank_in_species(loop), trim(temp), trim(elnes_label(edge_species(loop)))
      end if
    end do

    ! Now we know how many edges we have we can write them to a file

    ! Open the output file
    core_unit = io_file_unit()
    open (unit=core_unit, action='write', file=trim(seedname)//'_core_edge.dat')

    ! Write into the output file
    write (core_unit, *) '#*********************************************'
    write (core_unit, *) '#            Core loss function               '
    write (core_unit, *) '#*********************************************'
    write (core_unit, *) '#'
    if (core_LAI_broadening) then
      if (LAI_gaussian) write (core_unit, *) '# Gaussian broadening: FWHM', LAI_gaussian_width
      if (LAI_lorentzian) then
        write (core_unit, *) '# Lorentzian broadening included'
        write (core_unit, *) '# Lorentzian scale ', LAI_lorentzian_scale
        write (core_unit, *) '# Lorentzian offset ', LAI_lorentzian_offset
        write (core_unit, *) '# Lorentzian width ', LAI_lorentzian_width
      end if
    end if
    write (core_unit, *) '#'

    do loop = 1, num_edge
      write (core_unit, *) '# ', trim(edge_name(loop))

      dos_temp = 0.0_dp; dos_temp2 = 0.0_dp

      ! Have had to reallocate this in order to do the core_chemical_shift below
      if (set_efermi_zero) then
        E_shift = E - efermi
      else
        E_shift = E
      end if

      if (nspins == 1) then
        do loop2 = 1, edge_num_am(loop)
          dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos(:, 1, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
          if (core_LAI_broadening) then
            dos_temp2(:, 1) = dos_temp2(:, 1) + weighted_dos_broadened(:, 1, edge_list(loop, loop2))&
            &/real(edge_num_am(loop), dp)
          end if
        end do
      else
        do loop2 = 1, edge_num_am(loop)
          dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos(:, 1, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
          dos_temp(:, 2) = dos_temp(:, 2) + weighted_dos(:, 2, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
          dos_temp(:, 3) = dos_temp(:, 3) + dos_temp(:, 1) + dos_temp(:, 2)
          if (core_LAI_broadening) then
            dos_temp2(:, 1) = dos_temp2(:, 1) + weighted_dos_broadened(:, 1, edge_list(loop, loop2))/&
            &real(edge_num_am(loop), dp)
            dos_temp2(:, 2) = dos_temp2(:, 2) + weighted_dos_broadened(:, 2, edge_list(loop, loop2))/&
            &real(edge_num_am(loop), dp)
            dos_temp2(:, 3) = dos_temp2(:, 3) + dos_temp2(:, 1) + dos_temp2(:, 2)
          end if
        end do
      end if

      !Originally written to calculate the first nonzero term in the edge
      !elnes_edge = 0.0_dp

      !do N = 1, dos_nbins !doing this because we want to find the last 0.0000 value before the start of the peak edge
      !  if (dos_temp(N, 1) > 0.0_dp) then
      !    elnes_edge = E_shift(N)
      !    exit
      !  end if
      !end do

      !write (core_unit, *) elnes_edge !test to write out elnes_edge
      !write (core_unit, *) cbm_energy ! test to see if cbm calculated
      !write (core_unit, *) vbm_energy! test to see if cbm calculated
      ! Applies mizoguchi chemical shift if added to dos
      if (core_chemical_shift /= -1.0_dp) then
        E_shift = E + core_chemical_shift - cbm_energy
      end if

      do N = 1, dos_nbins
        if (nspins == 1) then
          if (core_LAI_broadening) then
            write (core_unit, *) E_shift(N), dos_temp(N, 1), dos_temp2(N, 1)
          else
            write (core_unit, *) E_shift(N), dos_temp(N, 1)
          end if
        else
          if (core_LAI_broadening) then
            write (core_unit, '(7(E21.13,2x))') E_shift(N), dos_temp(N, 1), dos_temp(N, 2),&
            & dos_temp(N, 3), dos_temp2(N, 1), dos_temp2(N, 2), dos_temp2(N, 3)
          else
            write (core_unit, '(4(E21.13,2x))') E_shift(N), dos_temp(N, 1), dos_temp(N, 2), dos_temp(N, 3)
          end if
        end if
      end do

      write (core_unit, *) ''
    end do

    close (core_unit)

    if (num_sites == 1) then ! if only one site we write out plot script files

      if (trim(output_format) == "xmgrace") then

        core_unit = io_file_unit()
        open (unit=core_unit, file=trim(seedname)//'_'//'core_edge'//'.agr', iostat=ierr)
        if (ierr .ne. 0) call io_error(" ERROR: Cannot open xmgrace batch file in core: write_core_xmgrace")

        min_x = minval(E_shift)
        max_x = maxval(E_shift)

        min_y = minval(weighted_dos)
        max_y = maxval(weighted_dos)

        ! For aesthetic reasons we make the axis range 1% larger than the data range
        range = abs(max_y - min_y)
        max_y = max_y + 0.01_dp*range
        min_y = 0.0_dp!min_y-0.01_dp*range

        call xmgu_setup(core_unit)
        call xmgu_legend(core_unit)
        call xmgu_title(core_unit, min_x, max_x, min_y, max_y, 'Core-loss Spectrum')

        call xmgu_axis(core_unit, "y", 'Units')
        call xmgu_axis(core_unit, "x", 'Energy (eV)')

        do loop = 1, num_edge
          dos_temp = 0.0_dp; dos_temp2 = 0.0_dp
          do loop2 = 1, edge_num_am(loop)
            if (nspins == 1) then
              dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos(:, 1, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
            else
              dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos(:, 1, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
              dos_temp(:, 2) = dos_temp(:, 2) + weighted_dos(:, 2, edge_list(loop, loop2))/real(edge_num_am(loop), dp)
              dos_temp(:, 3) = dos_temp(:, 1) + dos_temp(:, 1) + dos_temp(:, 2)
            end if
          end do

          call xmgu_data_header(core_unit, loop, loop, trim(edge_name(loop)))
          call xmgu_data(core_unit, loop, E_shift(:), dos_temp(:, 1)) ! only 1 element at the moment RJN 28Aug14
        end do

        if (core_LAI_broadening) then

          core_unit = io_file_unit()
          open (unit=core_unit, file=trim(seedname)//'_'//'core_edge_broad'//'.agr', iostat=ierr)
          if (ierr .ne. 0) call io_error(" ERROR: Cannot open xmgrace batch file in core: write_core_xmgrace")

          min_x = minval(E_shift)
          max_x = maxval(E_shift)

          min_y = minval(weighted_dos)
          max_y = maxval(weighted_dos)

          ! For aesthetic reasons we make the axis range 1% larger than the data range
          range = abs(max_y - min_y)
          max_y = max_y + 0.01_dp*range
          min_y = 0.0_dp!min_y-0.01_dp*range

          call xmgu_setup(core_unit)
          call xmgu_legend(core_unit)
          call xmgu_title(core_unit, min_x, max_x, min_y, max_y, 'Core-loss Spectrum')

          call xmgu_axis(core_unit, "y", 'Units')
          call xmgu_axis(core_unit, "x", 'Energy (eV)')

          do loop = 1, num_edge
            dos_temp = 0.0_dp; dos_temp2 = 0.0_dp
            do loop2 = 1, edge_num_am(loop)
              if (nspins == 1) then
                dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos_broadened(:, 1, edge_list(loop, loop2))/&
                &real(edge_num_am(loop), dp)
              else
                dos_temp(:, 1) = dos_temp(:, 1) + weighted_dos_broadened(:, 1, edge_list(loop, loop2))/&
                &real(edge_num_am(loop), dp)
                dos_temp(:, 2) = dos_temp(:, 2) + weighted_dos_broadened(:, 2, edge_list(loop, loop2))/&
                &real(edge_num_am(loop), dp)
                dos_temp(:, 3) = dos_temp(:, 1) + dos_temp(:, 1) + dos_temp(:, 2)
              end if
            end do

            call xmgu_data_header(core_unit, loop, loop, trim(edge_name(loop)))
            call xmgu_data(core_unit, loop, E_shift(:), dos_temp(:, 1))   ! Only 1 part for now RJN 28Aug14
          end do

        end if

      end if

    end if

    deallocate (E_shift, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating E_shift in write_core')

  end subroutine write_core

  subroutine core_gaussian
    !**************************************************************
    ! This subroutine adds in instrumental (Gaussian) broadening

    use od_constants, only: pi, dp, bohr2ang
    use od_parameters, only: LAI_gaussian_width, dos_nbins, LAI_lorentzian, LAI_lorentzian_scale
    use od_dos_utils, only: E
    use od_electronic, only: nspins, elnes_mwab, band_energy

    integer :: N, N_spin, N_energy, N_energy2
    real(kind=dp) :: G_width, g, dE
    real(kind=dp), allocatable, dimension(:, :, :) :: weighted_dos_temp

    G_width = LAI_gaussian_width          ! FWHM of Gaussian
    dE = E(2) - E(1)
    allocate (weighted_dos_temp(size(weighted_dos, 1), size(weighted_dos, 2), size(weighted_dos, 3)))
    weighted_dos_temp = 0.0_dp

    if (LAI_lorentzian .or. (LAI_lorentzian_scale .gt. 0.00001_dp)) then
      weighted_dos_temp = weighted_dos_broadened           ! In case we've already done Lorentzian broadening
      weighted_dos_broadened = 0.0_dp
    else
      weighted_dos_temp = weighted_dos
    end if

    do N = 1, elnes_mwab%norbitals         ! Loop over orbitals
      do N_spin = 1, nspins               ! Loop over spins
        do N_energy = 1, dos_nbins       ! Loop over energy
          do N_energy2 = 1, dos_nbins   ! Turn each energy value into a function
            g = (((4.0_dp*log(2.0_dp))/pi)**(0.5_dp))*(1/G_width)*exp(-4.0_dp*(log(2.0_dp))* &
                                                                      (((E(N_energy2) - E(N_energy))/G_width)**2.0_dp))  ! Gaussian
            weighted_dos_broadened(N_energy2, N_spin, N) = weighted_dos_broadened(N_energy2, N_spin, N) &
                                                           + (g*weighted_dos_temp(N_energy, N_spin, N)*dE)
          end do
        end do                        ! End loop over energy
      end do                           ! End loop over spins
    end do                              ! End loop over orbitals

  end subroutine core_gaussian

  subroutine core_lorentzian
    !**************************************************************
    ! This subroutine adds in life-time (Lorentzian) broadening

    use od_constants, only: pi, dp
    use od_parameters, only: LAI_lorentzian_width, LAI_lorentzian_scale, LAI_lorentzian_offset, &
      LAI_gaussian_width, dos_nbins, LAI_gaussian, adaptive, linear, fixed
    use od_dos_utils, only: E
    use od_electronic, only: nspins, elnes_mwab, efermi
    use od_dos_utils, only: efermi_fixed, efermi_adaptive, efermi_linear

    integer :: N, N_spin, N_energy, N_energy2
    real(kind=dp) :: L_width, l, dE

    dE = E(2) - E(1)

    do N = 1, elnes_mwab%norbitals         ! Loop over orbitals
      do N_spin = 1, nspins               ! Loop over spins
        do N_energy = 1, dos_nbins       ! Loop over energy
          if (E(N_energy) .ge. (LAI_lorentzian_offset + efermi)) then
            L_width = 0.5_dp*(LAI_lorentzian_width & ! HWHW of Lorentzian
                              + ((E(N_energy) - efermi - LAI_lorentzian_offset)*LAI_lorentzian_scale))
          else
            L_width = 0.5_dp*LAI_lorentzian_width
          end if
          if ((L_width*pi) .lt. dE) then  ! to get rid of spikes caused by L_width too small
            L_width = dE/pi
          end if
          do N_energy2 = 1, dos_nbins ! Turn each energy value into a function
            l = weighted_dos(N_energy, N_spin, N)*L_width/(pi*(((E(N_energy2) - E(N_energy))**2) + (L_width**2)))  ! Lorentzian
            weighted_dos_broadened(N_energy2, N_spin, N) = weighted_dos_broadened(N_energy2, N_spin, N) + (l*dE)
          end do
          !  end if
        end do                        ! End look over energy
      end do                           ! End loop over spins
    end do                              ! End loop over orbitals

  end subroutine core_lorentzian

!!$
!!$  !===============================================================================
!!$  subroutine write_core_gnuplot(label,E,column1,column2,column3)
!!$    !===============================================================================
!!$    use od_io,         only : io_file_unit,io_error,seedname
!!$    implicit none
!!$
!!$    type(graph_labels),intent(in) :: label
!!$
!!$    real(dp),  intent(in) :: E(:)
!!$    real(dp),  intent(in)  :: column1(:)
!!$    real(dp),  optional, intent(in) :: column2(:)
!!$    real(dp),  optional, intent(in) :: column3(:)
!!$
!!$    integer :: gnu_unit,ierr
!!$
!!$    gnu_unit=io_file_unit()
!!$    open(unit=gnu_unit,file=trim(seedname)//'_'//trim(label%name)//'.gnu',iostat=ierr)
!!$    if(ierr.ne.0) call io_error(" ERROR: Cannot open gnuplot batch file in optics: write_optics_gnupot")
!!$
!!$    gnu_unit = io_file_unit()
!!$    open(unit=gnu_unit,action='write',file=trim(seedname)//'_'//trim(label%name)//'.gnu')
!!$    write(gnu_unit,*) 'set xlabel ','"'//trim(label%x_label)//'"'
!!$    write(gnu_unit,*) 'set ylabel ','"'//trim(label%y_label)//'"'
!!$    write(gnu_unit,*) 'set title ','"'//trim(label%title)//'"'
!!$    if(present(column3)) then
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:3 t ','"'//trim(label%legend_b)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:4 t ','"'//trim(label%legend_c)//'"',' w l'
!!$    elseif(present(column2)) then
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l, \'
!!$       write(gnu_unit,*) '       "'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:3 t ','"'//trim(label%legend_b)//'"',' w l'
!!$    else
!!$       write(gnu_unit,*) 'plot ','"'//trim(seedname)//'_'//trim(label%name)//'.dat'//'"',' u 1:2 t ','"'//trim(label%legend_a)//'"',' w l'
!!$    endif
!!$    close(gnu_unit)
!!$
!!$  end subroutine write_optics_gnuplot
!!$

end module od_core
