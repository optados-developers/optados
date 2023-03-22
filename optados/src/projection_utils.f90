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
! MODULE od_projection_utils
! This module contains all of the routines to read projections of DOS and dispersion
! curves from CASTEP pdos_bin files, as well as parsing the user-requested projectors.
! It provides most calculated arrays as globals, namely:
!
!    - projection_array: array of projections for all projectors returned by CASTEP.
!    - matrix_weights: weights of desired projectors by (projector, band, kpoint, spin).
!    - proj_{symbol,am,sites}: array of projections split by OptaDOS shortcuts.
!
! as well as a few global variables, `num_proj` and `shortcut`.
!-------------------------------------------------------------------------------
module od_projection_utils

  use od_constants, only: dp
  use od_parameters, only: iprint
  use od_io, only: stdout, stderr

  implicit none

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  real(kind=dp), allocatable, public, dimension(:, :, :, :) :: matrix_weights
  integer, public, allocatable :: projection_array(:, :, :, :)
  integer, public :: num_proj

  integer, public, parameter :: max_am = 4                   ! s,p,d,f  hard coded!

  ! Data derived from the info in the pdos_weights file
  character(len=8), public, allocatable :: proj_symbol(:)   ! symbols
  integer, public, allocatable :: proj_am(:, :)              ! angular mtm (num_species,max_am)
  integer, public, allocatable :: proj_sites(:)             ! number of each species
  logical, public :: shortcut
  !-------------------------------------------------------------------------!

  character(len=1), parameter :: c_atomsep = ":"
  character(len=1), parameter :: c_labelsep = ":"
  character(len=10), parameter :: c_digit = "0123456789"
  character(len=1), parameter :: c_range = "-"
  character(len=1), parameter :: c_projsep = ","
  character(len=1), parameter :: c_delimiter = '"'
  character(len=4), parameter :: c_punc = " ,-:"

  private

  public :: projection_merge, projection_get_string, projection_analyse_orbitals

contains

  !===============================================================================
  subroutine projection_merge
    !===============================================================================
    ! This subroutine accumulates the weights of desired projectors into the
    ! matrix_weights array.
    !===============================================================================
    use od_electronic, only: pdos_orbital, pdos_weights, pdos_mwab, nspins
    use od_cell, only: num_kpoints_on_node
    use od_comms, only: my_node_id
    use od_io, only: io_error, stdout
    implicit none

    integer :: N, N_spin, n_eigen, nproj, orb, ierr

    allocate (matrix_weights(num_proj, pdos_mwab%nbands, num_kpoints_on_node(my_node_id), nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_merge - allocation of matrix_weights failed')
    matrix_weights = 0.0_dp

    do N = 1, num_kpoints_on_node(my_node_id) ! Loop over kpoints
      do N_spin = 1, nspins                  ! Loop over spins
        do n_eigen = 1, pdos_mwab%nbands    ! Loop over unoccupied states
          do nproj = 1, num_proj
            do orb = 1, pdos_mwab%norbitals
              if (projection_array(pdos_orbital%species_no(orb), pdos_orbital%rank_in_species(orb) &
                                   , pdos_orbital%am_channel(orb) + 1, nproj) == 1) then
                matrix_weights(nproj, n_eigen, N, N_spin) = matrix_weights(nproj, n_eigen, N, N_spin) + &
                                                            pdos_weights(orb, n_eigen, N, N_spin)
              end if
            end do
          end do
        end do
      end do
    end do

    return

  end subroutine projection_merge

  !===============================================================================
  subroutine projection_get_string
    !! Read the PDOS string and pouplate the projection_array.
    !! AJM pre Feb 2022: This is a mindbendingly horrific exercise in book-keeping

    use od_parameters, only: projectors_string
    use od_cell, only: num_species, atoms_species_num, atoms_label
    use od_io, only: maxlen, io_error
    use od_algorithms, only: channel_to_am
    implicit none
    character(len=maxlen) :: ctemp, ctemp2, ctemp3

    integer :: loop4, loop3, loop2, ierr
    logical :: pdos_sum
    logical :: long_atom_name_exits, delimiter_exists
    integer   :: loop, pos, loop_l, loop_a, loop_p
    integer   ::  species_count, species_proj

    integer, allocatable :: pdos_temp(:, :, :, :)

    !Check for any short cuts
    shortcut = .false.

    ctemp = projectors_string
    if (index(ctemp, 'species_ang') > 0) then
      num_proj = 0
      do loop = 1, num_species
        num_proj = num_proj + count(proj_am(loop, :) == 1)
      end do
      allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')

      projection_array = 0
      loop_p = 1
      do loop = 1, num_species
        do loop_l = 1, max_am
          if (proj_am(loop, loop_l) == 0) cycle
          projection_array(loop, :, loop_l, loop_p) = 1
          loop_p = loop_p + 1
        end do
      end do
      shortcut = .true.
    elseif (index(ctemp, 'species') > 0) then
      num_proj = num_species
      allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')

      projection_array = 0
      do loop = 1, num_species
        projection_array(loop, :, :, loop) = 1
      end do
      shortcut = .true.
    elseif (index(ctemp, 'sites') > 0) then
      num_proj = 0
      do loop = 1, num_species
        num_proj = num_proj + proj_sites(loop)
      end do
      allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')

      projection_array = 0
      loop_p = 1
      do loop = 1, num_species
        do loop_a = 1, proj_sites(loop)
          projection_array(loop, loop_a, :, loop_p) = 1
          loop_p = loop_p + 1
        end do
      end do
      shortcut = .true.
    elseif (index(ctemp, 'angular') > 0) then
      num_proj = max_am
      allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')
      projection_array = 0
      loop_p = 0
      do loop = 1, num_proj
        projection_array(:, :, loop, loop) = 1
      end do
      shortcut = .true.
    end if

    if (.not. shortcut) then

      !look for sum
      ctemp = projectors_string
      pdos_sum = .false.
      if (index(ctemp, 'sum:') == 1) then
        pdos_sum = .true.
        ctemp = ctemp(5:)
      end if
      ! take 1st part of string

      ctemp2 = ctemp
      species_count = 0; num_proj = 0
      do
        delimiter_exists = .false.
        long_atom_name_exits = .false.
        !! look for each species section
        !! and pasrse to find number of projections
        call projection_find_atom(ctemp2, ctemp3, delimiter_exists, pos)
        call projection_analyse_atom(ctemp3, long_atom_name_exits, species_proj)
        if (long_atom_name_exits .and. .not. delimiter_exists) then
          call io_error('Error: projection_get_string - something odd reading the PDOS string. BUG (2)')
        elseif (.not. long_atom_name_exits .and. delimiter_exists) then
          if (iprint > 2) write (stdout, *) " projection_get_string: delimiter but &
          &no atom label. Not invalid. Just unnecessary."
        end if
        num_proj = num_proj + species_proj
        species_count = species_count + 1
        if (pos == 0 .or. pos > len(trim(ctemp2))) exit ! There's nothing left to read
        ctemp2 = ctemp2(pos + 1:)
      end do

      ! now allocate the correct sized array
      allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
      if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')
      projection_array = 0

      ctemp2 = ctemp
      do loop = 1, species_count
        !! loop for each species section and
        !! then fill projection_array
        call projection_find_atom(ctemp2, ctemp3, delimiter_exists, pos)
        call projection_analyse_atom(ctemp3, long_atom_name_exits)
        if (long_atom_name_exits .and. .not. delimiter_exists) then
          call io_error('Error: projection_get_string - something odd reading the PDOS string. BUG (2)')
        elseif (.not. long_atom_name_exits .and. delimiter_exists) then
          if (iprint > 2) write (stdout, *) " projection_get_string: delimiter but &
          &no atom label. Not invalid. Just unnecessary."
        end if
        ctemp2 = ctemp2(pos + 1:)
      end do

      if (pdos_sum) then
        allocate (pdos_temp(num_species, maxval(atoms_species_num), max_am, 1), stat=ierr)
        if (ierr /= 0) call io_error('Error: projection_get_string - allocation of pdos_temp failed')

        pdos_temp = 0
        do loop4 = 1, num_proj
          do loop3 = 1, max_am
            do loop2 = 1, maxval(atoms_species_num)
              do loop = 1, num_species
                if (projection_array(loop, loop2, loop3, loop4) == 1) then
                  pdos_temp(loop, loop2, loop3, 1) = 1
                end if
              end do
            end do
          end do
        end do
        deallocate (projection_array, stat=ierr)
        if (ierr /= 0) call io_error('Error: projection_get_string - deallocation of projection_array failed')
        num_proj = 1
        allocate (projection_array(num_species, maxval(atoms_species_num), max_am, num_proj), stat=ierr)
        if (ierr /= 0) call io_error('Error: projection_get_string - allocation of projection_array failed')
        projection_array = 0
        projection_array = pdos_temp
      end if
    end if

    if (iprint > 2) then
      write (stdout, *) "+--------------------------------------------------------------------------+"
      write (stdout, *) "|  projection_array( num_species, atoms_in_species, max_am, num_proj)      |"
      write (stdout, *) "|       Atoms         Species     Atoms in                Array            |"
      write (stdout, *) "| Proj  label          Number     Species      am         value            |"
      do loop4 = 1, num_proj
        write (stdout, *) "+--------------------------------------------------------------------------+"
        do loop3 = 1, max_am
          do loop2 = 1, maxval(atoms_species_num)
            do loop = 1, num_species
              write (stdout, '(1x,a1,i4,6x,a10,4x,i4,6x,i4,8x,1a,8x,i4,15x,a1)') "|", loop4, atoms_label(loop),&
              & loop, loop2, channel_to_am(loop3), projection_array(loop, loop2, loop3, loop4), "|"
            end do
          end do
        end do
      end do
      write (stdout, *) "+--------------------------------------------------------------------------+"
    end if

    return

  end subroutine projection_get_string

  !===============================================================================
  subroutine projection_find_atom(cstring_in, catom_out, atom_label, atomsep_position)
    !! Looks for the first atom in cstring_in and returns it in catom_out
    !! returns the position of the seperator afer the atom or zero if there is none
    !! after it.
    !! Returns atom label to indicate that this is not just a one or two character
    !! atom name.

    use od_io, only: maxlen, io_error, stdout

    character(len=maxlen), intent(in) :: cstring_in
    character(len=maxlen), intent(out) :: catom_out
    integer, intent(out) :: atomsep_position
    logical, intent(out) :: atom_label
    integer :: delimiter_position, delimiter_position_end

    if (iprint > 2) then
      write (stdout, *) "+-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*+"
      write (stdout, *) "|                 projection_untils :>  projection_find_atom               |"
      write (stdout, *) "+--------------------------------------------------------------------------+"
      write (stdout, '(1x,a1,a30,1x,a20,23x,a1)') "|", "  String in :", trim(cstring_in), '|'
    end if

    if (len(trim(cstring_in)) == 0) call io_error('projection_find_atom: Zero lenght string input. This is a BUG!')

    atom_label = .false.

    atomsep_position = index(cstring_in, c_atomsep)
    delimiter_position = index(cstring_in, c_delimiter)

    ! If the position is zero there is no seperator, so the whole string is the
    ! atom
    if ((atomsep_position == 0) .and. (delimiter_position == 0)) then
      ! There's one atom and no delimiting
      if (iprint > 2) write (stdout, '(1x,a1,a30,44x,a1)') "|", " one atom, no delimiter.", "|"
      catom_out = cstring_in
    elseif (delimiter_position .ne. 0) then
      ! Look for the second delimiter
      delimiter_position_end = index(cstring_in(delimiter_position + 1:), c_delimiter)
      if (delimiter_position_end == 0) call io_error('projection_find_atom: Error &
          & finding closing inverted comma')
      ! The seperator is now after the end delimiter
      atomsep_position = index(cstring_in(delimiter_position_end:), c_atomsep)
      if (atomsep_position == 0) then
        catom_out = cstring_in
      else
        ! The atom seperator is after the delimiter
        atomsep_position = atomsep_position + delimiter_position_end - 1
        catom_out = cstring_in(:atomsep_position - 1)
      end if
      atom_label = .true.
      if (iprint > 2) write (stdout, '(1x,a1,a30,1x,i20,23x,a1)') "|", " atom seperator position :", atomsep_position, "|"
      if (iprint > 2) write (stdout, '(1x,a1,a30,1x,i10,1x,i10,22x,a1)') "|", "  delimiters at position :",&
      & delimiter_position, delimiter_position_end, "|"
    elseif ((atomsep_position .ne. 0) .and. (delimiter_position == 0)) then
      ! If there aren't any delimiters than we just jump to the end of the atom seperator
      if (iprint > 2) write (stdout, *) " no delimiters"
      catom_out = cstring_in(:atomsep_position - 1)
    else
      call io_error('projection_find_atom: Error parsing string. This is a BUG!')
    end if

    catom_out = trim(catom_out)
    if (iprint > 2) then
      if (atom_label) write (stdout, '(1x,a1,a30,1x,a20,23x,a1)') "|", "Atom label present  :", "TRUE", '|'
      write (stdout, '(1x,a1,a30,1x,a20,23x,a1)') "|", "  String out :", trim(catom_out), '|'
      write (stdout, *) "+-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*+"
      write (stdout, *)
    end if
    return

  end subroutine projection_find_atom

  !===============================================================================
  subroutine projection_analyse_atom(ctemp, atom_label, species_proj)
    !! Receive a string about an indivudial species in ctemp and either count the
    !! number of projectors needed (if species_proj is present) or actually
    !! populate the module varaible projection_array(:, :, :, :) (if species_proj)
    !! is absent.
    !! This subroutine is called twice, the first time it counts the number of
    !! projectors necessary such that projection_array(:, :, :, :) can be
    !! allocated. The second time it populates projection_array(:, :, :, :)
    !! AJM Feb 2022
    !! AJM (Before Feb 2022): This is a mindbendingly horrific exercise in book-keeping

    use od_cell, only: num_species, atoms_species_num, atoms_label, atoms_symbol
    use od_io, only: maxlen, io_error
    implicit none

    character(len=maxlen), intent(inout) :: ctemp
    !! Provides a snippet of a PDOS string that contains info about a
    !! particular species or projector e.g. "C:exi"1(s,p)

    logical, intent(out) :: atom_label
    !! Whether the atom we're looking at has a label or not e.g. C:exi where
    !! exi is the label

    integer, optional, intent(out) :: species_proj
    !! Count the number of projectors we need for this snippet of the PDOS string

    integer, allocatable :: pdos_atoms(:), pdos_ang(:)

    integer, save :: offset = 0
    !! Keep tabs on how may projectors we have already allocated. We call this
    !! function multiple times as we build up the module variable
    !! projection_array(:, :, :, :)

    character(len=3)  :: c_symbol = '   '
    !! The atomic symbol of the species we have in the snippet

    character(len=maxlen) ::  c_am
    !! The angular momentum information stripped off the input snippet

    character(len=maxlen) ::  catom_label
    !! If present, the label after the atom's species in the input snippet.
    !! i.e. the "exi" in C:exi.

    character(len=maxlen) ::  cspecies
    !! A temporary variable contining the species name whilst we're stripping
    !! out the delimiters.

    character(len=maxlen) :: m_string
    !! String containing an individual momentum channel label, i.e. s,p,f

    character(len=maxlen) :: dummy, label, ctemp2
    !! Dummy labels

    character(len=5)  :: c_num1, c_num2

    integer :: species, num_sites, num_am

    integer :: i_digit, ispecies

    integer :: pos_l, pos_r, ia, iz, idiff, ic1, ic2

    logical :: am_sum, site_sum

    integer   :: num1, num2, i_punc, pos3, loop_l, loop_a, loop_p, loop_j
    integer   :: counter, loop_r, range_size, ierr, label_position

    integer :: delimiter_position_start, delimiter_position_end
    !! Start and end of where the quotes are, if any.

    logical :: lcount
    !! If true we are just counting projectors, if false we have already allocated
    !! the PDOS arrays are looking to populate them.

    logical :: delimiter_exists
    !! If the atom label has quotes aroud it

    allocate (pdos_atoms(maxval(atoms_species_num)), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_analyse_atom - allocation of pdos_atoms failed')
    allocate (pdos_ang(max_am), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_analyse_atom - allocation of pdos_ang failed')

    pdos_atoms = 0; pdos_ang = 0
    site_sum = .false.
    am_sum = .false.
    atom_label = .false.
    delimiter_exists = .false.
    lcount = .false.

    if (present(species_proj)) lcount = .true.
    !! We are just in projector counting mode. Otherwise we're in
    !! projection_array(:, :, :, :) filling mode.

    ctemp = trim(ctemp)

    if (iprint > 2) then
      write (stdout, *) "  +-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+"
      write (stdout, *) "  |             projection_untils :>  projection_analyse_atom            |"
      write (stdout, *) "  +----------------------------------------------------------------------+"
      write (stdout, '(3x,a1,a30,1x,a20,19x,a1)') "|", "  String in :", trim(ctemp), '|'
    end if

    !! Parse for Ang Mtm string eg (s,p). The angular momentum information is
    !! is stripped off into c_am and ctemp goes forward without it.
    c_am = ''
    pos_l = index(ctemp, '(')
    if (pos_l > 0) then
      ! There are brackets
      pos_r = index(ctemp, ')')
      if (pos_r == 0) call io_error('projection_analyse_atom: found ( but no )')
      if (pos_r <= pos_l) call io_error('projection_analyse_atom: found ) before (')
      c_am = ctemp(pos_l + 1:pos_r - 1)
      ctemp = ctemp(:pos_l - 1)
    else ! implicit sum over AM
      am_sum = .true.
    end if

    if (iprint > 2) write (stdout, '(3x,a1,a30,1x,a20,19x,a1)') "|", "  String with AM removed :", trim(ctemp), '|'

    !! Parse for inverted commas around the species name. If they are present
    !! we strip them off and send ctemp forward without them.
    !! We set delimiter_exist = .true. so that we know that they were there.
    !! Note: if a label is present, we send that on as well, s.t. "C:exi" ==> C:exi
    !! the routines below will cope with the colon (:) it just had to be delimited
    !! in this subroutine's parent routine.
    !! Note that we could delimit without a label, e.g. "C", we do this s.t. later
    !! we only match C and not with C:exi.
    delimiter_position_start = index(ctemp, c_delimiter)
    if (delimiter_position_start .ne. 0) then
      !! We have a delimiter. Need to find the other end.
      delimiter_position_end = index(ctemp(delimiter_position_start + 1:), c_delimiter) + delimiter_position_start
      if (iprint > 2) write (stdout, '(3x,a1,a30,1x,i10,1x,i10,18x,a1)') "|", "  delimiters at position :",&
      & delimiter_position_start, delimiter_position_end, "|"
      if (delimiter_position_end == 0) call io_error('projection_analyse_atom: Opening " found but no closing "')
      !! Strip the atom symbol from the symbol:label
      catom_label = ctemp(delimiter_position_start + 1:delimiter_position_end)
      !! Remove the symbols for the delimiter -- but keep the rest intact in ctemp
      ctemp = ctemp(delimiter_position_start + 1:delimiter_position_end - 1) &
      &//trim(ctemp(delimiter_position_end + 1:))
      delimiter_exists = .true.
    else
      !! We do not have a delimiter
      catom_label = ''
      delimiter_exists = .false.
    end if

    !! Parse for an atom label, e.g. C:exi. ctemp goes forward with the
    !! label stripped off. If a label exists cspecies and catom_label are populated and
    !! atom_label is set to true.
    label_position = index(catom_label, c_labelsep)
    if (label_position == 0) then
      ! There isn't a label seperator
      atom_label = .false. ! Say it explicitly
      catom_label = ''
      if (iprint > 2) write (stdout, '(3x,a1,a30,1x,a4,35x,a1)') "|", "  Label found :", "None", '|'
      ctemp = trim(ctemp)
    elseif (label_position > 1) then
      ! There is a seperator
      atom_label = .true.
      catom_label = ctemp(label_position + 1:delimiter_position_end - 2)
      if (iprint > 2) write (stdout, '(3x,a1,a30,1x,a20,19x,a1)') "|", "  Label found :", trim(catom_label), '|'
      cspecies = trim(ctemp(1:label_position - 1))
      if (iprint > 2) write (stdout, '(3x,a1,a30,1x,a20,19x,a1)') "|", "  Species found :", trim(cspecies), '|'
      ctemp = trim(cspecies)//ctemp(delimiter_position_end - 1:) !Subtract 1 extra becasue we've stripped the ""
    else ! Error
      write (stderr, *) 'projection_analyse_atom: cannot understand atom &
      &label, ', ctemp
      call io_error('projection_analyse_atom: found an atom label but could &
      &not work out a valid syntax for it')
    end if

    !! Now we parse the atomic symbol placing it into c_symbol. Whatever is left
    !! is sent forward in ctemp. In this case that can only be the number in species.
    ia = ichar('a')
    iz = ichar('z')
    idiff = ichar('Z') - ichar('z')

    ic1 = ichar(ctemp(1:1))
    if (ic1 < ia .or. ic1 > iz) call io_error('projection_analyse_atom: problem reading atomic symbol in pdos string')
    ic2 = ichar(ctemp(2:2))
    if (ic2 >= ia .and. ic1 <= iz) then
      c_symbol(1:1) = char(ic1 + idiff)
      c_symbol(2:2) = ctemp(2:2)
      ctemp = ctemp(3:)
    else
      c_symbol(1:1) = char(ic1 + idiff)
      c_symbol(2:2) = ''
      ctemp = ctemp(2:)
    end if

    !! I don't understand this code. It seems convoluted. I've reproduced it below
    !! to cope with atom labels.
    !! The code finds the species number for the species in the PDOS snippet.
    !! I think it is for checking whether the species in the string matches any
    !! species in the pdos_bin
!      species = 0
!      do loop_j = 1, num_species
!        if (adjustl(c_symbol) == adjustl(proj_symbol(loop_j))) then
!          species = loop_j
!        end if
!      end do
    species = 0
    do loop_j = 1, num_species
      if (.not. atom_label) then
        if (trim(c_symbol) == adjustl(proj_symbol(loop_j))) then
          species = loop_j
        end if
      elseif (atom_label) then
        if (trim(c_symbol)//":"//trim(catom_label) == adjustl(proj_symbol(loop_j))) then
          species = loop_j
        end if
      end if
    end do

    if (iprint > 2) then
      write (stdout, '(3x,a1,a30,1x,i5,34x,a1)') "|", "     Species index number  :", species, "|"
      write (stdout, '(3x,a1,a30,1x,a5,34x,a1)') "|", "     Atom number in species :", trim(ctemp), "|"
    end if

    !! Check for not matching a species
    if (species == 0) then
      if (.not. atom_label) then
        write (stderr, *) "PDOS string atom ", trim(c_symbol), " matches none of ", adjustl(proj_symbol(:))
      elseif (atom_label) then
        write (stderr, *) "PDOS string atom ", trim(c_symbol)//":"//trim(catom_label), " matches none of ", adjustl(proj_symbol(:))
      endif
      call io_error('projection_analyse_atom: Failed to match atomic symbol in pdos string')
    endif

    !! Parse if there's a particular species site number. ie. He3
    counter = 0
    dummy = adjustl(ctemp)
    if (len_trim(dummy) > 0) then
      dummy = adjustl(dummy)
      do
        i_punc = scan(dummy, c_punc)
        if (i_punc == 0) call io_error('projection_analyse_atom: error looking for atom numbers')
        c_num1 = dummy(1:i_punc - 1)
        read (c_num1, *, err=101, end=101) num1
        dummy = adjustl(dummy(i_punc:))
        !look for range
        if (scan(dummy, c_range) == 1) then
          i_digit = scan(dummy, c_digit)
          dummy = adjustl(dummy(i_digit:))
          i_punc = scan(dummy, c_punc)
          c_num2 = dummy(1:i_punc - 1)
          read (c_num2, *, err=101, end=101) num2
          write (stdout, *) "dummy=", dummy
          dummy = adjustl(dummy(i_punc:))
          range_size = abs(num2 - num1) + 1
          do loop_r = 1, range_size
            counter = counter + 1
            if (min(num1, num2) + loop_r - 1 > proj_sites(species)) &
               call io_error('projection_analyse_atom: Atom number given in pdos string &
       &is greater than number of atoms for given species')
            pdos_atoms(min(num1, num2) + loop_r - 1) = 1
          end do
        else
          counter = counter + 1
          if (num1 > proj_sites(species)) &
             call io_error('projection_analyse_atom: Atom number given in pdos string &
        &is greater than number of atoms for given species')
          pdos_atoms(num1) = 1
        end if

        if (scan(dummy, c_projsep) == 1) dummy = adjustl(dummy(2:))
        if (scan(dummy, c_range) == 1) &
             & call io_error('projection_analyse_atom: Error parsing atoms numbers - incorrect range')
        if (index(dummy, ' ') == 1) exit
      end do
    else
      site_sum = .true.
    end if

    !! Parse the angular momentum part of the input snippet that we found earlier
    !! and placed in c_am
    counter = 0
    dummy = adjustl(c_am)
    if (iprint > 2) write (stdout, '(3x,a1,a30,1x,a20,19x,a1)') "|", " Ang. mom. channels found :", trim(dummy), "|"
    if (len_trim(dummy) > 0) then
      do
        pos3 = index(dummy, ',')
        if (pos3 == 0) then
          ctemp2 = dummy
        else
          ctemp2 = dummy(:pos3 - 1)
        end if
        read (ctemp2(1:), *, err=106, end=106) m_string
        select case (trim(adjustl(m_string)))
        case ('s')
          pdos_ang(1) = 1
        case ('p')
          pdos_ang(2) = 1
        case ('d')
          pdos_ang(3) = 1
        case ('f')
          pdos_ang(4) = 1
        case default
          call io_error('projection_analyse_atom: Problem reading l state ')
        end select
        if (pos3 == 0) exit
        dummy = dummy(pos3 + 1:)
      end do
    else
      am_sum = .true.
    end if

    !! If we have particularly labelled sites (e.g. He3) then we need a projector
    !! for each site. Otherwise we just need a projector for the atom
    if (site_sum) then
      num_sites = 1
    else
      num_sites = count(pdos_atoms == 1)
    end if

    !! If we have particularly labelled A.M. channels (e.g. C(s,p)) then  we
    !! need a projector for each site. Otherwise we just need a projector for
    !! the atom.
    if (am_sum) then
      num_am = 1
    else
      num_am = count(pdos_ang == 1)
    end if

    !! If all we're supposed to be doing is counting projectors. Then count them
    !! and leave. Else try to populate the projection_array
    if (lcount) then
      species_proj = num_am*num_sites
    elseif (.not. lcount) then
      loop_p = 1 + offset
      if (site_sum .and. am_sum) then
        projection_array(species, :, :, loop_p) = 1
        if (.not. delimiter_exists .and. catom_label == '') then
          !! Its not in quotes and its not
          !! specially tagged. Hence C must match any C:exi present
          do ispecies = 1, num_species
            if (ispecies == species) cycle !self interaction!
            if (atoms_symbol(ispecies) == c_symbol) then
              projection_array(ispecies, :, :, loop_p) = 1
            end if
          end do
        end if
      elseif (site_sum .and. .not. am_sum) then
        do loop_l = 1, max_am
          if (pdos_ang(loop_l) == 0) cycle
          projection_array(species, :, loop_l, loop_p) = 1
          if (.not. delimiter_exists .and. catom_label == '') then
            !! Its not in quotes and its not
            !! specially tagged. Hence C must match C:exi
            do ispecies = 1, num_species
              if (ispecies == species) cycle !self interaction!
              if (atoms_symbol(ispecies) == c_symbol) then
                projection_array(ispecies, :, loop_l, loop_p) = 1
              end if
            end do
          end if
          loop_p = loop_p + 1
        end do
      elseif (.not. site_sum .and. am_sum) then
        do loop_a = 1, proj_sites(species)
          if (pdos_atoms(loop_a) == 0) cycle
          projection_array(species, loop_a, :, loop_p) = 1
          loop_p = loop_p + 1
        end do
      else
        do loop_l = 1, max_am
          if (pdos_ang(loop_l) == 0) cycle
          do loop_a = 1, proj_sites(species)
            if (pdos_atoms(loop_a) == 0) cycle
            projection_array(species, loop_a, loop_l, loop_p) = 1
            loop_p = loop_p + 1
          end do
        end do
      end if
      offset = loop_p - 1

    end if

    if (iprint > 2) then
      write (stdout, *) "  +-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=+"
      write (stdout, *)
    end if

    return

101 call io_error('projection_analyse_atom Error parsing keyword ')
106 call io_error('projection_analyse_atom: Problem reading l state into string ')

  end subroutine projection_analyse_atom

  subroutine projection_analyse_orbitals
    !! Find out what projectors we have read in.
    use od_electronic, only: pdos_orbital, pdos_mwab
    use od_cell, only: atoms_symbol, num_species, atoms_label
    use od_constants, only: periodic_table_name
    use od_io, only: io_error
    implicit none

    integer :: loop, loop2, counter, ierr
    !! Only loop counters and temp variables

    if (iprint > 2) then
      write (stdout, *) "+--------------------------------------------------------------------------+"
      write (stdout, *) "|                 projection_untils :>  projection_analyse_orbitals        |"
      write (stdout, *) "+--------------------------------------------------------------------------+"
    end if

    if (maxval(pdos_orbital%species_no(:)) > num_species) &
      call io_error('Error: projection_analyse_substring - more species in pdos file than in cell file')

    allocate (proj_sites(maxval(pdos_orbital%species_no(:))), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_analyse_substring - allocation of proj_sites failed')
    allocate (proj_am(maxval(pdos_orbital%species_no(:)), max_am), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_analyse_substring - allocation of proj_am failed')
    allocate (proj_symbol(maxval(pdos_orbital%species_no(:))), stat=ierr)
    if (ierr /= 0) call io_error('Error: projection_analyse_substring - allocation of proj_symbol failed')
    proj_sites = 0; proj_am = 0
    do loop = 1, pdos_mwab%norbitals
      if (pdos_orbital%rank_in_species(loop) > proj_sites(pdos_orbital%species_no(loop))) &
        proj_sites(pdos_orbital%species_no(loop)) = pdos_orbital%rank_in_species(loop)
      if (pdos_orbital%rank_in_species(loop) == 1) &
        proj_am(pdos_orbital%species_no(loop), pdos_orbital%am_channel(loop) + 1) = 1
    end do

    ! Now need to figure out symbols for each species

! AJM: This looks like it doesn't put things in CASTEP atom order if labels are present
! TBH I don't understand why it was doing it this way -- will probably learn the hard way! :)
!    counter = 1
!    do loop2 = 1, 109
!      do loop = 1, num_species
!        if (atoms_symbol(loop) == periodic_table_name(loop2)) then
!          proj_symbol(counter) = periodic_table_name(loop2)
!          counter = counter + 1
!          !check atom count here
!        end if
!      end do
!    end do
!! This matches the atomic_symbols from the cell file to the proj_symbol of the PDOS
!! routines
    counter = 0
    do loop = 1, num_species
      do loop2 = 1, 109
        if (atoms_symbol(loop) == periodic_table_name(loop2)) then
          proj_symbol(loop) = atoms_label(loop)
          counter = counter + 1
        end if
      end do
    end do

    if (counter .ne. num_species) call io_error('Error: projection_analyse_substring &
    &- mismatch in atoms from .pdos_bin and out.cell')

    if (iprint > 2) then
      write (stdout, '(1x,a1,1x,10x,1x,a14,3x,a15,3x,a14,13x,a1)') "|", "atoms_label(:)", &
      &"atoms_symbol(:)", "proj_symbol(:)", "|"
      do loop = 1, num_species
        write (stdout, '(1x,a1,1x,i5,1x,a14,3x,a15,3x,a14,18x,a1)') "|", loop, trim(atoms_label(loop)),&
        &trim(atoms_symbol(loop)), trim(proj_symbol(loop)), "|"
      end do

      do loop = 1, num_species
        do loop2 = loop + 1, num_species
          if (proj_symbol(loop) == proj_symbol(loop2)) Then
            write (stdout, '(1x,a1,1x,a24,3x,a3,3x,a3,3x,i5,3x,i5,21x,a1)') "|", "duplicate species found:", &
            &proj_symbol(loop), proj_symbol(loop2), loop, loop2, "|"
          end if
        end do
      end do
      write (stdout, *) "+--------------------------------------------------------------------------+"
    end if
  end subroutine projection_analyse_orbitals

end module od_projection_utils
