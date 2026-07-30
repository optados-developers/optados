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
! Module: Cell                                                            !
! For dealing with the real and reciprocal space cell                     !
!-------------------------------------------------------------------------!
! Modules used:  constants                                                !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below                                                         !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris (So far)                            11/10/2010 !
!=========================================================================!
module od_cell
  use od_constants, only: dp
  use od_io, only: maxlen
  implicit none

  private
  !-------------------------------------------------------------------------!
  ! R E A L   S P A C E   V A R I A B L E S
  real(kind=dp), public, save :: real_lattice(1:3, 1:3)
  real(kind=dp), public, save :: recip_lattice(1:3, 1:3)
  real(kind=dp), public, save :: cell_volume
  !-------------------------------------------------------------------------!

  !-------------------------------------------------------------------------!
  ! R E C I P R O C A L   S P A C E   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: kpoint_r(:, :)
  real(kind=dp), allocatable, public, save :: kpoint_r_cart(:, :)
  real(kind=dp), allocatable, public, save :: kpoint_weight(:)
  integer, allocatable, public, save       :: num_kpoints_on_node(:)

  integer, public, save :: nkpoints
  integer, public, save :: kpoint_grid_dim(3)

  !-------------------------------------------------------------------------!
  ! Symmetry Operations
  integer, public, save :: num_crystal_symmetry_operations
  real(kind=dp), allocatable, public, save :: crystal_symmetry_disps(:, :)
  real(kind=dp), allocatable, public, save :: crystal_symmetry_operations(:, :, :)

  ! Atom sites
  real(kind=dp), allocatable, public, save :: atoms_pos_frac(:, :, :)
  real(kind=dp), allocatable, public, save :: atoms_pos_cart(:, :, :)
  integer, allocatable, public, save :: atoms_species_num(:)
  character(len=maxlen), allocatable, public, save :: atoms_label(:)
  character(len=2), allocatable, public, save :: atoms_symbol(:)
  integer, public, save :: num_atoms
  integer, public, save :: num_species
  character(len=maxlen), allocatable, public, save  :: atoms_label_tmp(:)

  ! Added for photoemission
  real(kind=dp), allocatable, public, save :: atoms_pos_cart_photo(:, :)

  !-------------------------------------------------------------------------!
  ! G L O B A L L Y   A V A I L A B L E   F U N C T I O N S
  public :: cell_find_MP_grid
  public :: cell_calc_lattice
  public :: cell_report_parameters
  public :: cell_get_atoms
  public :: cell_read_cell
  public :: cell_get_symmetry
  public :: cell_dist
  ! Added for photoemission
  public :: cell_get_real_lattice
  public :: cell_calc_kpoint_r_cart
  !-------------------------------------------------------------------------!

contains

  !=========================================================================!
  subroutine cell_find_MP_grid(kpoints, num_kpts, kpoint_grid_dim, kpoint_offset)
    ! WARNING the kpoint_offset is only +/- the true kpoint offset. We need to
    ! do some more work to find out its sign
    ! A J Morris 29th September 2011
    !=========================================================================!
    use od_io, only: io_error, stdout
    implicit none
    integer, intent(out):: kpoint_grid_dim(1:3)
    integer, intent(in) :: num_kpts
    real(kind=dp), intent(in)     :: kpoints(1:3, 1:num_kpts)
    real(kind=dp), intent(out), optional :: kpoint_offset(1:3)
    real(kind=dp)            :: kpoint_TR(1:3, 1:num_kpts*2)
    real(kind=dp)            :: unique_kpoints(1:3, num_kpts*2)
    integer                  :: nunique_kpoints, iunique_kpoints
    integer                  :: ikpt, idim, jkpt, i
    real(kind=dp)    :: subtraction_tol, min_img, min_img2
    real(kind=dp)    :: min_img3, image, min_img_tol

    integer :: iprint = 1 ! We can't use the global iprint as it's higher that this
    ! in the module heirarchy.

    ! Before time reversal symmetry we could any combintaion of kpoints
    !--------------------X-------------
    !                  0.25

    ! We apply TR to make we have a complete set (and add periodic images)
    !------X------X------X------X------
    !    -0.75  -0.25   0.25  0.75
    !
    ! In this case the answer is easy. The minimim image is 0.5. Hence the grid is
    ! 1/0.5 = 2x MP.

    ! However the fun comes when shifts are applied. In this case applying TR makes
    ! things more complicated. Consider the same grid with a +ve 0.05 shift.
    !
    ! Now CASTEP had to give us more kpoints
    !------------------X------X------
    !                0.20   0.70
    !
    ! We couldn't know that an offset had been applied so we appy TR and get
    !   -0.70     -0.20  0.30    0.80
    !----X-X-----X-X----X-X----X-X----
    !  -0.80  -0.30  0.20   0.70
    !
    ! Hence now the 1st minimum image is 0.10, double the offset
    !           the 2nd    "      "   is 0.40
    !           the 3rd    "      "   is 0.50  the reciprocal of the MP grid

    ! One can now image the case where the offset is 3/(4n) where n is the MP grid number.
    ! in this case 3/8=0.375
    !                 (shifted -0.75)  (shifted -0.25)
    !------------------------X-------------X-------------------------
    !                     -0.375        0.125
    !
    ! After TR
    !                           -0.125        0.375
    !------------------------X-----X------X-----X--------------------
    !                     -0.375        0.125
    !
    ! So these are all equally spaced so algorithm sees this like the first example, easy, 0.25. Hence
    ! grid is 1/0.25 = 4x MP. Wrong!
    !
    ! This is a known bug -- and without having the symmetry ops, we can't build the shifted 2xMP grid to
    ! show that it is not a 4x MP.  AJM 15/9/2014

    ! When two numbers are the same
    subtraction_tol = epsilon(subtraction_tol)
    !write(*,*) "subtraction_tol=", subtraction_tol
    ! When one number is larger than another one
    min_img_tol = 0.000001_dp

    if (iprint > 3) then
      write (stdout, *)
      write (stdout, *) "+----------------------------------------------------------------------------+"
      write (stdout, *) "                              MP grid finder "
      write (stdout, *)
      write (stdout, *) " Kpoints found:", num_kpts
    end if

    ! Add time reversal
    kpoint_TR(1:3, 1:num_kpts) = kpoints(1:3, 1:num_kpts)
    kpoint_TR(1:3, num_kpts + 1:num_kpts*2) = -kpoints(1:3, 1:num_kpts)

    ! Act on each dimension independently
    do idim = 1, 3
      ! Fold all kpoints between (-0.5,0.5]
      do ikpt = 1, num_kpts*2
        kpoint_TR(idim, ikpt) = kpoint_TR(idim, ikpt) - floor(kpoint_TR(idim, ikpt) + 0.5_dp)
      end do
    end do

    unique_kpoints = 0.0_dp

    over_dim: do idim = 1, 3
      ! Make unique
      if (iprint > 3) write (stdout, *) " ----- Dimension ", idim, " -----"

      nunique_kpoints = 1
      unique_kpoints(idim, 1) = kpoint_TR(idim, 1)

      over_kpts: do ikpt = 2, 2*num_kpts
        do iunique_kpoints = 1, nunique_kpoints
          if (abs(unique_kpoints(idim, iunique_kpoints) - kpoint_TR(idim, ikpt)) .le. subtraction_tol) then
            !We've seen this before
            cycle over_kpts
          end if
        end do
        ! If we ended up here then, this is new
        nunique_kpoints = nunique_kpoints + 1
        if (iprint > 3) write (stdout, *) " ikpt= ", ikpt, "nunique_kpoints= ", nunique_kpoints
        unique_kpoints(idim, nunique_kpoints) = kpoint_TR(idim, ikpt)
      end do over_kpts

      if (iprint > 3) write (stdout, *) " Number of unique kpoints:", nunique_kpoints

      ! write(*,*) "------------------------ KPOINTS IN+TR+FOLDING+UNIQUE -----------"
      ! do i=1,nunique_kpoints
      !    write(*,*) i, idim, unique_kpoints(idim,i)
      ! enddo
      ! write(*,*) "-----------------------------------------------------------"
      if (iprint > 3) write (stdout, *) " Looking for special cases for dimension:", idim

      ! Look at special cases. These are largely necessary because the general finder needs to get to
      ! second nearest neighbour before it can work correctly. Except in the case of
      if (nunique_kpoints == 1) then
        ! There is only 1 k-point
        ! The shift is it's position (either 0 or 0.5)
        if (present(kpoint_offset)) kpoint_offset(idim) = unique_kpoints(idim, 1)
        kpoint_grid_dim(idim) = 1
        cycle over_dim
      elseif (nunique_kpoints == 2) then
        min_img = abs(unique_kpoints(idim, 1) - unique_kpoints(idim, 2))
        if (abs(min_img - 0.5_dp) .le. min_img_tol) then
          ! If the 1stMI is 0.5 then there are 2 k-points
          ! The shift is it's position (either 0 or 0.25)
          if (present(kpoint_offset)) kpoint_offset(idim) = unique_kpoints(idim, 1)
          kpoint_grid_dim(idim) = 2
          cycle over_dim
        else
          ! If the 1stMI .ne.0.5 then there is 1 kpoint
          ! It's offset is min_img/2
          if (present(kpoint_offset)) kpoint_offset(idim) = min_img/2
          kpoint_grid_dim(idim) = 1
          cycle over_dim
        end if
      elseif (nunique_kpoints == 3) then
        ! This is the case of a MP3 grid with a point at Gamma
        !! AJM COMMENTED OUT AS A TEST AGAINST
        !! HEXAGONAL CELLS
        !      if(present(kpoint_offset)) kpoint_offset(idim)=0.0_dp
        !      kpoint_grid_dim(idim)=3
        !      cycle over_dim
        continue
      end if

      if (iprint > 3) write (stdout, *) " Left the special-case kpoint finder for dimension:", idim
      if (iprint > 3) write (stdout, *) " Now trying to use the general solver..."

      ! Get 1st, 2nd and 3rd minimum images
      ! Get 1st min image
      min_img = huge(min_img)
      do ikpt = 1, nunique_kpoints
        do jkpt = ikpt + 1, nunique_kpoints
          image = abs(unique_kpoints(idim, ikpt) - unique_kpoints(idim, jkpt))
          if (image < min_img) min_img = image
        end do
      end do
      if (abs(min_img - huge(min_img)) .le. subtraction_tol) &
           & call io_error('cell_find_MP_grid: Failed to find a 1st min image')

      ! Get 2nd min image
      min_img2 = huge(min_img2)
      do ikpt = 1, nunique_kpoints
        do jkpt = ikpt, nunique_kpoints
          image = abs(unique_kpoints(idim, ikpt) - unique_kpoints(idim, jkpt))
          if ((image < min_img2) .and. image > min_img + min_img_tol) min_img2 = image
        end do
      end do
      if (abs(min_img2 - huge(min_img2)) .le. subtraction_tol) &
           & call io_error('cell_find_MP_grid: Failed to find a 2nd min image')

      ! Get 3rd min image
      min_img3 = huge(min_img3)
      do ikpt = 1, nunique_kpoints
        do jkpt = ikpt, nunique_kpoints
          image = abs(unique_kpoints(idim, ikpt) - unique_kpoints(idim, jkpt))
          if ((image < min_img3) .and. image > min_img2 + min_img_tol) min_img3 = image
        end do
      end do
      if (abs(min_img3 - huge(min_img3)) .le. subtraction_tol) then
        ! It's possible in CASTEP 6 and 7 to have 3 kpoints only for a 5 MP grid
        ! the specieal case solver then got confused.
        if (abs(2.0_dp*min_img - min_img2) < min_img_tol) then
          ! If 1stMI==2ndMI then 1/1stMP is the grid density and we're still ok. Carry on.
          ! we'll catch this scenario further down the code.
          if (iprint > 3) write (stdout, *) " Couldn't find a third minimum image, but since we don't think"
          if (iprint > 3) write (stdout, *) " there's an offset, it looks like a", 1.0_dp/min_img, " kpoint grid"
          continue
        else
          call io_error('cell_find_MP_grid: Failed to find a 3rd min image')
        end if
      end if

      if (iprint > 3) then
        write (stdout, *) " Min images for dimension:", idim
        write (stdout, *) " 1st:", min_img
        write (stdout, *) " 2nd:", min_img2
        write (stdout, *) " 3rd:", min_img3
        write (stdout, *) " 1st^-1:", 1.0_dp/min_img
      end if

      if (abs(2.0_dp*min_img - min_img2) < min_img_tol) then
        ! If 1stMI==2ndMI then 1/1stMP is the grid density
        if (present(kpoint_offset)) kpoint_offset(idim) = 0.0_dp
        kpoint_grid_dim(idim) = int(1.0_dp/min_img)
        ! WARNING could also have a shifted grid with a perfect shift 3/(4n)
        ! this would be a known bug
      else
        ! If 1stMI.ne.2ndMI then 1/3rdMP is the grid density
        ! and 1stMI/2 is the shift
        if (present(kpoint_offset)) kpoint_offset(idim) = min_img/2.0_dp
        kpoint_grid_dim(idim) = int(1.0_dp/min_img3)
      end if

    end do over_dim

    if (iprint > 3) write (stdout, *) " Conclusion = kpoint_grid_dim: ", kpoint_grid_dim
    ! if(present(kpoint_offset))  write(*,*) "kpoint_offset= ",  kpoint_offset

    if (iprint > 3) &
         & write (stdout, *) "+----------------------------------------------------------------------------+"
  end subroutine cell_find_MP_grid
  !=========================================================================!

  !=========================================================================!
  subroutine cell_get_symmetry
    !=========================================================================!
    ! Read in the cell symmetries                                             !
    !-------------------------------------------------------------------------!
    ! Arguments: kpoints - an array of kpoints                                !
    !            num_kpts - size of the kpoint array                          !
    !-------------------------------------------------------------------------!
    ! Returns: kpint_grid_dim - the number of kpoints in each dimension       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: None                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:  None                                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! Described below                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: None
    !--------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------!
    ! JRY, April 2011                                                         !
    !=========================================================================!
    use od_comms, only: on_root, comms_bcast
    use od_io, only: filename_len, io_file_unit, seedname, io_error, stdout
    implicit none
    integer :: ierr, sym_file
    logical :: exists
    character(filename_len)     :: sym_filename

    !check if we already have the symmetries
    if (allocated(crystal_symmetry_operations)) return

    sym_file = io_file_unit()
    sym_filename = trim(seedname)//".sym"
    if (on_root) inquire (file=sym_filename, exist=exists)
    call comms_bcast(exists, 1)

    if (.not. exists) then
      if (on_root) then
        write (stdout, '(1x,78a)') '!--------------------------------- WARNING ----------------------------------!'
        write (stdout, '(1x,78a)') '!                Symmetry Operations file (.sym) not found                   !'
        write (stdout, '(1x,78a)') '!                       Proceeding without symmetry                          !'
        write (stdout, '(1x,78a)') '!----------------------------------------------------------------------------!'
        write (stdout, '(1x,a78)') '|                                                                            |'
      end if
      num_crystal_symmetry_operations = 0
      return
    end if

    if (on_root) then
      open (unit=sym_file, file=sym_filename, form='unformatted', err=100, status='old')

      read (sym_file) num_crystal_symmetry_operations
      if (num_crystal_symmetry_operations > 0) then
        allocate (crystal_symmetry_operations(3, 3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_get_symmetry")
        allocate (crystal_symmetry_disps(3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_get_symmetry")
        read (sym_file) crystal_symmetry_operations
        read (sym_file) crystal_symmetry_disps
      end if
    end if

    call comms_bcast(num_crystal_symmetry_operations, 1)
    if (num_crystal_symmetry_operations > 0) then
      if (.not. on_root) then
        allocate (crystal_symmetry_operations(3, 3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_get_symmetry")
        allocate (crystal_symmetry_disps(3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_get_symmetry")
      end if
      call comms_bcast(crystal_symmetry_operations(1, 1, 1), 9*num_crystal_symmetry_operations)
      call comms_bcast(crystal_symmetry_disps(1, 1), 3*num_crystal_symmetry_operations)
    end if

    return

100 call io_error('Error: Problem opening sym file in cell_get_symmetry')

  end subroutine cell_get_symmetry
  !=========================================================================!

  !=========================================================================!
  subroutine cell_read_cell
    !=========================================================================!
    use od_constants, only: bohr2ang
    use od_io, only: io_file_unit, io_error, seedname, maxlen
    use od_algorithms, only: utility_cart_to_frac, utility_frac_to_cart, utility_lowercase

    implicit none

    real(kind=dp), allocatable     :: atoms_pos_frac_tmp(:, :)
    real(kind=dp), allocatable    :: atoms_pos_cart_tmp(:, :)
    character(len=20) :: keyword
    integer           :: in, in1, in2, ins, ine, loop, i, line_e, line_s, counter, tot_num_lines
    integer           :: loop2, max_sites, ierr, ic, num_lines, line_counter, in_unit
    logical           :: found_e, found_s, frac
    character(len=maxlen) :: dummy
    character(len=maxlen), allocatable :: ctemp(:)
    !character(len=maxlen), allocatable :: atoms_label_tmp(:)
    logical           :: lconvert

    character(len=maxlen), allocatable :: in_data(:)

    ! read in the cell file

    ! count the lines
    in_unit = io_file_unit()
    open (in_unit, file=trim(seedname)//'-out.cell', form='formatted', status='old', err=101)

    num_lines = 0; tot_num_lines = 0
    do
      read (in_unit, '(a)', iostat=ierr, err=200, end=210) dummy
      dummy = adjustl(dummy)
      tot_num_lines = tot_num_lines + 1
      if (.not. dummy(1:1) == '!' .and. .not. dummy(1:1) == '#') then
        if (len(trim(dummy)) > 0) num_lines = num_lines + 1
      end if

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'-out.cell')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'-out.cell')
210 continue
    rewind (in_unit)

    ! now read in for real - ignoring comments
    allocate (in_data(num_lines), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating in_data in cell_get_atoms')

    line_counter = 0
    do loop = 1, tot_num_lines
      read (in_unit, '(a)', iostat=ierr, err=200) dummy
      dummy = utility_lowercase(dummy)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '!' .or. dummy(1:1) == '#') cycle
      if (len(trim(dummy)) == 0) cycle
      line_counter = line_counter + 1
      in1 = index(dummy, '!')
      in2 = index(dummy, '#')
      if (in1 == 0 .and. in2 == 0) in_data(line_counter) = dummy
      if (in1 == 0 .and. in2 > 0) in_data(line_counter) = dummy(:in2 - 1)
      if (in2 == 0 .and. in1 > 0) in_data(line_counter) = dummy(:in1 - 1)
      if (in2 > 0 .and. in1 > 0) in_data(line_counter) = dummy(:min(in1, in2) - 1)
    end do

    close (in_unit)

    ! let's look for the atoms block (remember everything is lower case)
    keyword = 'positions'

    found_s = .false.
    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), '%block')
      if (in == 0 .or. in > 1) cycle
      if (index(in_data(loop), 'frac') > 0) then
        frac = .true.
      elseif (index(in_data(loop), 'abs') > 0) then
        frac = .false.
      else
        cycle
      end if
      line_s = loop
      if (found_s) then
        call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
      end if
      found_s = .true.
    end do

    if (frac) then
      keyword = 'positions_frac'
    else
      keyword = 'positions_abs'
    end if

    found_e = .false.
    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), '%endblock')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
      end if
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found %block'//trim(keyword)//' but no %endblock'//trim(keyword)//' in cell file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: %endblock'//trim(keyword)//' comes before %block'//trim(keyword)//' in input file')
    end if

    ! now we know where the atoms block is

    lconvert = .false.
    dummy = in_data(line_s + 1)
    if (index(dummy, 'ang') .ne. 0) then
      lconvert = .false.
      line_s = line_s + 1
    elseif (index(dummy, 'bohr') .ne. 0) then
      lconvert = .true.
      line_s = line_s + 1
    end if

    num_atoms = line_e - 1 - (line_s + 1) + 1
    allocate (atoms_pos_frac_tmp(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac_tmp in cell_get_atoms')
    allocate (atoms_pos_cart_tmp(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart_tmp in cell_get_atoms')
    allocate (ctemp(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ctemp in cell_get_atoms')
    allocate (atoms_label_tmp(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label_tmp in cell_get_atoms')

    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (frac) then
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_frac_tmp(i, counter), i=1, 3)
      else
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_cart_tmp(i, counter), i=1, 3)
      end if
    end do

    if (lconvert) then
      atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr2ang
    end if

    call cell_get_real_lattice
    if (frac) then
      do loop = 1, num_atoms
        call utility_frac_to_cart(atoms_pos_frac_tmp(:, loop), atoms_pos_cart_tmp(:, loop), real_lattice)
      end do
    else
      do loop = 1, num_atoms
        call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), atoms_pos_frac_tmp(:, loop), recip_lattice)
      end do
    end if

    ! Now we sort the data into the proper structures
    num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          num_species = num_species + 1
          ctemp(num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms_species_num(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in cell_get_atoms')
    allocate (atoms_label(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in cell_get_atoms')
    allocate (atoms_symbol(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in cell_get_atoms')
    atoms_species_num(:) = 0

    do loop = 1, num_species
      atoms_label(loop) = ctemp(loop)
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms_species_num(loop) = atoms_species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms_species_num)
    allocate (atoms_pos_frac(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in cell_get_atoms')
    allocate (atoms_pos_cart(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in cell_get_atoms')
    allocate (atoms_pos_cart_photo(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart_photo in cell_get_atoms')

    ! Making a copy to use in the photo.f90 subroutine "analyse_geometry"
    atoms_pos_cart_photo = atoms_pos_cart_tmp

    do loop = 1, num_species
      counter = 0
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms_pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms_pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, num_species
      atoms_symbol(loop) (1:2) = atoms_label(loop) (1:2)
      ic = ichar(atoms_symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms_symbol(loop) (2:2) = ' '
    end do

    ! let's look for the symmetry_ops block
    keyword = 'symmetry_ops'
    num_crystal_symmetry_operations = 0
    found_s = .false.
    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), '%block')
      if (in == 0 .or. in > 1) cycle
      line_s = loop
      if (found_s) then
        call io_error('Error: Found %block'//trim(keyword)//' more than once in out.cell file')
      end if
      found_s = .true.
    end do

    if (found_s) then

      found_e = .false.
      do loop = 1, num_lines
        ine = index(in_data(loop), trim(keyword))
        if (ine == 0) cycle
        in = index(in_data(loop), '%endblock')
        if (in == 0 .or. in > 1) cycle
        line_e = loop
        if (found_e) then
          call io_error('Error: Found %block'//trim(keyword)//' more than once in out.cell file')
        end if
        found_e = .true.
      end do

      if (.not. found_e) then
        call io_error('Error: Found %block'//trim(keyword)//' but no %endblock'//trim(keyword)//' in out.cell file')
      end if

      if (line_e <= line_s) then
        call io_error('Error: %endblock'//trim(keyword)//' comes before %block'//trim(keyword)//' in input file')
      end if

      ! now we know where the block is
      num_crystal_symmetry_operations = (line_e - line_s - 1)/4
      if ((4*num_crystal_symmetry_operations) /= (line_e - line_s - 1)) &
        call io_error('Error: Something wrong with symmetry_ops block in -out.cell file')
      if (num_crystal_symmetry_operations > 0) then
        allocate (crystal_symmetry_operations(3, 3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_read_cell")
        allocate (crystal_symmetry_disps(3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_read_cell")
        counter = 0
        do loop = line_s + 1, line_e - 1, 4
          dummy = in_data(loop)
          counter = counter + 1
          dummy = in_data(loop)
          read (dummy, *, err=240, end=240) crystal_symmetry_operations(1, 1, counter), &
            crystal_symmetry_operations(2, 1, counter), crystal_symmetry_operations(3, 1, counter)
          dummy = in_data(loop + 1)
          read (dummy, *, err=240, end=240) crystal_symmetry_operations(1, 2, counter), &
            crystal_symmetry_operations(2, 2, counter), crystal_symmetry_operations(3, 2, counter)
          dummy = in_data(loop + 2)
          read (dummy, *, err=240, end=240) crystal_symmetry_operations(1, 3, counter), &
            crystal_symmetry_operations(2, 3, counter), crystal_symmetry_operations(3, 3, counter)
          dummy = in_data(loop + 3)
          read (dummy, *, err=240, end=240) crystal_symmetry_disps(1, counter), &
            crystal_symmetry_disps(2, counter), crystal_symmetry_disps(3, counter)
        end do
      end if
    else
      call io_error('Error: Cannot find %block '//trim(keyword)//' in '//trim(seedname)//'-out.cell')
    end if

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine cell_read_cell
  !=========================================================================!

  !=========================================================================!
  subroutine cell_get_atoms
    !=========================================================================!
    use od_constants, only: bohr2ang
    use od_io, only: io_file_unit, io_error, seedname, maxlen
    use od_algorithms, only: utility_cart_to_frac, utility_frac_to_cart, utility_lowercase

    implicit none

    real(kind=dp), allocatable     :: atoms_pos_frac_tmp(:, :)
    real(kind=dp), allocatable    :: atoms_pos_cart_tmp(:, :)
    character(len=20) :: keyword
    integer           :: in, in1, in2, ins, ine, loop, i, line_e, line_s, counter, tot_num_lines
    integer           :: loop2, max_sites, ierr, ic, num_lines, line_counter, in_unit
    logical           :: found_e, found_s, frac
    character(len=maxlen) :: dummy
    character(len=maxlen), allocatable :: ctemp(:)
    character(len=maxlen), allocatable :: atoms_label_tmp(:)
    logical           :: lconvert

    character(len=maxlen), allocatable :: in_data(:)

    ! read in the cell file

    ! count the lines
    in_unit = io_file_unit()
    open (in_unit, file=trim(seedname)//'.cell', form='formatted', status='old', err=101)

    num_lines = 0; tot_num_lines = 0
    do
      read (in_unit, '(a)', iostat=ierr, err=200, end=210) dummy
      dummy = adjustl(dummy)
      tot_num_lines = tot_num_lines + 1
      if (.not. dummy(1:1) == '!' .and. .not. dummy(1:1) == '#') then
        if (len(trim(dummy)) > 0) num_lines = num_lines + 1
      end if

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.cell')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.cell')
210 continue
    rewind (in_unit)

    ! now read in for real - ignoring comments
    allocate (in_data(num_lines), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating in_data in cell_get_atoms')

    line_counter = 0
    do loop = 1, tot_num_lines
      read (in_unit, '(a)', iostat=ierr, err=200) dummy
      dummy = utility_lowercase(dummy)
      dummy = adjustl(dummy)
      if (dummy(1:1) == '!' .or. dummy(1:1) == '#') cycle
      if (len(trim(dummy)) == 0) cycle
      line_counter = line_counter + 1
      in1 = index(dummy, '!')
      in2 = index(dummy, '#')
      if (in1 == 0 .and. in2 == 0) in_data(line_counter) = dummy
      if (in1 == 0 .and. in2 > 0) in_data(line_counter) = dummy(:in2 - 1)
      if (in2 == 0 .and. in1 > 0) in_data(line_counter) = dummy(:in1 - 1)
      if (in2 > 0 .and. in1 > 0) in_data(line_counter) = dummy(:min(in1, in2) - 1)
    end do

    close (in_unit)

    ! let's look for the atoms block (remember everything is lower case)
    keyword = 'positions'

    found_s = .false.
    do loop = 1, num_lines
      ins = index(in_data(loop), trim(keyword))
      if (ins == 0) cycle
      in = index(in_data(loop), '%block')
      if (in == 0 .or. in > 1) cycle
      if (index(in_data(loop), 'frac') > 0) then
        frac = .true.
      elseif (index(in_data(loop), 'abs') > 0) then
        frac = .false.
      else
        cycle
      end if
      line_s = loop
      if (found_s) then
        call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
      end if
      found_s = .true.
    end do

    if (frac) then
      keyword = 'positions_frac'
    else
      keyword = 'positions_abs'
    end if

    found_e = .false.
    do loop = 1, num_lines
      ine = index(in_data(loop), trim(keyword))
      if (ine == 0) cycle
      in = index(in_data(loop), '%endblock')
      if (in == 0 .or. in > 1) cycle
      line_e = loop
      if (found_e) then
        call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
      end if
      found_e = .true.
    end do

    if (.not. found_e) then
      call io_error('Error: Found %block'//trim(keyword)//' but no %endblock'//trim(keyword)//' in cell file')
    end if

    if (line_e <= line_s) then
      call io_error('Error: %endblock'//trim(keyword)//' comes before %block'//trim(keyword)//' in input file')
    end if

    ! now we know where the atoms block is

    lconvert = .false.
    dummy = in_data(line_s + 1)
    if (index(dummy, 'ang') .ne. 0) then
      lconvert = .false.
      line_s = line_s + 1
    elseif (index(dummy, 'bohr') .ne. 0) then
      lconvert = .true.
      line_s = line_s + 1
    end if

    num_atoms = line_e - 1 - (line_s + 1) + 1
    allocate (atoms_pos_frac_tmp(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac_tmp in cell_get_atoms')
    allocate (atoms_pos_cart_tmp(3, num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart_tmp in cell_get_atoms')
    allocate (ctemp(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating ctemp in cell_get_atoms')
    allocate (atoms_label_tmp(num_atoms), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label_tmp in cell_get_atoms')

    counter = 0
    do loop = line_s + 1, line_e - 1
      dummy = in_data(loop)
      counter = counter + 1
      if (frac) then
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_frac_tmp(i, counter), i=1, 3)
      else
        read (dummy, *, err=240, end=240) atoms_label_tmp(counter), (atoms_pos_cart_tmp(i, counter), i=1, 3)
      end if
    end do

    if (lconvert) then
      atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr2ang
    end if

    call cell_get_real_lattice

    if (frac) then
      do loop = 1, num_atoms
        call utility_frac_to_cart(atoms_pos_frac_tmp(:, loop), atoms_pos_cart_tmp(:, loop), real_lattice)
      end do
    else
      do loop = 1, num_atoms
        call utility_cart_to_frac(atoms_pos_cart_tmp(:, loop), atoms_pos_frac_tmp(:, loop), recip_lattice)
      end do
    end if

    ! Now we sort the data into the proper structures
    num_species = 1
    ctemp(1) = atoms_label_tmp(1)
    do loop = 2, num_atoms
      do loop2 = 1, loop - 1
        if (trim(atoms_label_tmp(loop)) == trim(atoms_label_tmp(loop2))) exit
        if (loop2 == loop - 1) then
          num_species = num_species + 1
          ctemp(num_species) = atoms_label_tmp(loop)
        end if
      end do
    end do

    allocate (atoms_species_num(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_species_num in cell_get_atoms')
    allocate (atoms_label(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_label in cell_get_atoms')
    allocate (atoms_symbol(num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_symbol in cell_get_atoms')
    atoms_species_num(:) = 0

    do loop = 1, num_species
      atoms_label(loop) = ctemp(loop)
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          atoms_species_num(loop) = atoms_species_num(loop) + 1
        end if
      end do
    end do

    max_sites = maxval(atoms_species_num)
    allocate (atoms_pos_frac(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in cell_get_atoms')
    allocate (atoms_pos_cart(3, max_sites, num_species), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in cell_get_atoms')

    do loop = 1, num_species
      counter = 0
      do loop2 = 1, num_atoms
        if (trim(atoms_label(loop)) == trim(atoms_label_tmp(loop2))) then
          counter = counter + 1
          atoms_pos_frac(:, counter, loop) = atoms_pos_frac_tmp(:, loop2)
          atoms_pos_cart(:, counter, loop) = atoms_pos_cart_tmp(:, loop2)
        end if
      end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop = 1, num_species
      atoms_symbol(loop) (1:2) = atoms_label(loop) (1:2)
      ic = ichar(atoms_symbol(loop) (2:2))
      if ((ic .lt. ichar('a')) .or. (ic .gt. ichar('z'))) &
        atoms_symbol(loop) (2:2) = ' '
    end do

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine cell_get_atoms

  !=========================================================================!
  subroutine cell_calc_lattice
    !=========================================================================!
    ! Begin with a real lattice. Convert from bohr. Calculate a reciprocal
    ! lattice and the volume of the cell
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables used: real_lattice, recip_lattice, cell_volume
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by Andrew Morris from the LinDOS program             11/10/2010
    !=========================================================================
    use od_constants, only: pi, bohr2ang
    implicit none

    ! THESE ARE IN BOHR, DON'T GET TRIPPED UP AGAIN!
    real_lattice = real_lattice*bohr2ang

    call cell_get_real_lattice

    recip_lattice(1, 1) = real_lattice(2, 2)*real_lattice(3, 3) - &
                          real_lattice(3, 2)*real_lattice(2, 3)
    recip_lattice(2, 1) = real_lattice(2, 3)*real_lattice(3, 1) - &
                          real_lattice(3, 3)*real_lattice(2, 1)
    recip_lattice(3, 1) = real_lattice(2, 1)*real_lattice(3, 2) - &
                          real_lattice(3, 1)*real_lattice(2, 2)
    recip_lattice(1, 2) = real_lattice(3, 2)*real_lattice(1, 3) - &
                          real_lattice(1, 2)*real_lattice(3, 3)
    recip_lattice(2, 2) = real_lattice(3, 3)*real_lattice(1, 1) - &
                          real_lattice(1, 3)*real_lattice(3, 1)
    recip_lattice(3, 2) = real_lattice(3, 1)*real_lattice(1, 2) - &
                          real_lattice(1, 1)*real_lattice(3, 2)
    recip_lattice(1, 3) = real_lattice(1, 2)*real_lattice(2, 3) - &
                          real_lattice(2, 2)*real_lattice(1, 3)
    recip_lattice(2, 3) = real_lattice(1, 3)*real_lattice(2, 1) - &
                          real_lattice(2, 3)*real_lattice(1, 1)
    recip_lattice(3, 3) = real_lattice(1, 1)*real_lattice(2, 2) - &
                          real_lattice(2, 1)*real_lattice(1, 2)

    ! * Calculate cell volume
    cell_volume = real_lattice(1, 1)*recip_lattice(1, 1) + &
                  real_lattice(2, 1)*recip_lattice(1, 2) + &
                  real_lattice(3, 1)*recip_lattice(1, 3)

    if (cell_volume < 0.0_dp) then ! Left handed set
      cell_volume = -cell_volume
    end if

    ! Scale reciprocal lattice by 2*pi/volume
    recip_lattice(:, :) = recip_lattice(:, :)*pi*2.0_dp/cell_volume

  end subroutine cell_calc_lattice

  !=========================================================================!
  subroutine cell_report_parameters
    !=========================================================================
    ! Begin with a real lattice. Convert from bohr. Calculate a reciprocal
    ! lattice and the volume of the cell
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables used: real_lattice, recip_lattice, cell_volume
    !-------------------------------------------------------------------------
    ! Modules used:  See below
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------
    ! Necessary conditions: None
    !-------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------
    ! Written by J R Yates, modified A J Morris                     Dec 2010
    !=========================================================================
    use od_io, only: stdout

    implicit none

    integer :: i

    write (stdout, '(30x,a21)') 'Lattice Vectors (Ang)'
    write (stdout, 101) 'a_1', (real_lattice(1, I), i=1, 3)
    write (stdout, 101) 'a_2', (real_lattice(2, I), i=1, 3)
    write (stdout, 101) 'a_3', (real_lattice(3, I), i=1, 3)

    write (stdout, *)

    write (stdout, '(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    write (stdout, 101) 'b_1', (recip_lattice(1, I), i=1, 3)
    write (stdout, 101) 'b_2', (recip_lattice(2, I), i=1, 3)
    write (stdout, 101) 'b_3', (recip_lattice(3, I), i=1, 3)

    write (stdout, *)
    write (stdout, '(19x,a17,3x,f11.5)', advance='no') &
      'Unit Cell Volume:', cell_volume

    write (stdout, '(2x,a7)') '(Ang^3)'
    write (stdout, *)

    return
101 format(20x, a3, 2x, 3F11.6)

  end subroutine cell_report_parameters

  subroutine cell_dist
    use od_comms, only: comms_bcast, on_root
    use od_io, only: io_error
    implicit none

    integer :: max_sites, ierr

    call comms_bcast(real_lattice(1, 1), 9)
    call comms_bcast(recip_lattice(1, 1), 9)
    call comms_bcast(cell_volume, 1)

    !call comms_bcast(kpoint_r(:,:)
    !call comms_bcast(kpoint_r_cart(:,:)
    !call comms_bcast(kpoint_weight(:)

    call comms_bcast(nkpoints, 1)
    call comms_bcast(kpoint_grid_dim(1), 3)

    !-------------------------------------------------------------------------!

    call comms_bcast(num_atoms, 1)
    call comms_bcast(num_species, 1)
    if (num_atoms > 0) then
      if (.not. on_root) then
        allocate (atoms_species_num(num_species), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_species_num in cell_dist')
      end if
      call comms_bcast(atoms_species_num(1), num_species)
      max_sites = maxval(atoms_species_num)
      if (.not. on_root) then
        allocate (atoms_pos_frac(3, max_sites, num_species), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_pos_frac in cell_dist')
        allocate (atoms_pos_cart(3, max_sites, num_species), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_pos_cart in cell_dist')
        allocate (atoms_label(num_species), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_label in cell_dist')
        allocate (atoms_symbol(num_species), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_symbol in cell_dist')
        ! For Photoemission
        allocate (atoms_pos_cart_photo(3, num_atoms), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_pos_cart_photo in cell_dist')
        allocate (atoms_label_tmp(num_atoms), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating atoms_label_tmp in cell_dist')
      end if
      call comms_bcast(atoms_pos_frac(1, 1, 1), 3*num_species*max_sites)
      call comms_bcast(atoms_pos_cart(1, 1, 1), 3*num_species*max_sites)
      call comms_bcast(atoms_label(1), len(atoms_label(1))*num_species)
      call comms_bcast(atoms_symbol(1), len(atoms_symbol(1))*num_species)
      ! For Photoemission
      call comms_bcast(atoms_pos_cart_photo(1, 1), 3*num_atoms)
      call comms_bcast(atoms_label_tmp(1), maxlen*num_atoms)
    end if
    call comms_bcast(num_crystal_symmetry_operations, 1)
    if (num_crystal_symmetry_operations > 0) then
      if (.not. on_root) then
        allocate (crystal_symmetry_operations(3, 3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_dist")
        allocate (crystal_symmetry_disps(3, num_crystal_symmetry_operations), stat=ierr)
        if (ierr /= 0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_dist")
      end if
      call comms_bcast(crystal_symmetry_operations(1, 1, 1), 3*3*num_crystal_symmetry_operations)
      call comms_bcast(crystal_symmetry_disps(1, 1), 3*num_crystal_symmetry_operations)
    end if

  end subroutine cell_dist

  !=========================================================================!
  subroutine cell_get_real_lattice
    !=========================================================================
    ! This subroutine reads the lattice parameters from the bands file in
    ! order to have them stored when the frac_to_cart and cart_frac subroutines
    ! are called. Independently of having used the elec_read_band_energy
    ! subroutine before.

    use od_comms, only: on_root
    use od_io, only: io_file_unit, seedname, filename_len, stdout, io_time, &
      io_error
    use od_constants, only: bohr2ang

    integer :: band_unit
    character(filename_len) :: band_filename

    !Open the bands file
    band_unit = io_file_unit()
    band_filename = trim(seedname)//".bands"
!    print*,'band_filename=',band_filename

    ! Read the header from the bands file
    if (on_root) then
      open (unit=band_unit, file=band_filename, status="old", form='formatted')!,err=100)
!100    call io_error('Error: Problem opening bands file in cell_get_real_lattice')
      read (band_unit, *)
      read (band_unit, *)
      read (band_unit, *)
      read (band_unit, *)
      read (band_unit, *)
      read (band_unit, *)
      read (band_unit, *) real_lattice(:, 1)
      read (band_unit, *) real_lattice(:, 2)
      read (band_unit, *) real_lattice(:, 3)
    end if
    real_lattice = real_lattice*bohr2ang
    if (on_root) close (unit=band_unit)

  end subroutine cell_get_real_lattice

  !=========================================================================!
  subroutine cell_calc_kpoint_r_cart
    !=========================================================================
    ! This subroutine calculates the cartesian coordinates of the k points
    use od_algorithms, only: utility_reciprocal_frac_to_cart
    use od_comms, only: my_node_id
    use od_io, only: io_file_unit, seedname, filename_len, stdout, io_time, &
      io_error
!    use od_electronic, only : elec_read_band_energy

    integer :: i, ik, loop, ierr
    real(kind=dp), allocatable, dimension(:, :) :: kpoint_r_tmp
    real(kind=dp), allocatable, dimension(:, :) :: kpoint_r_cart_tmp

    allocate (kpoint_r_tmp(3, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating kpoint_r_tmp in&
&    cell_calc_kpoint_r_cart')
    allocate (kpoint_r_cart_tmp(3, num_kpoints_on_node(my_node_id)), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating kpoint_r_cart_tmp in&
&    cell_calc_kpoint_r_cart')

    kpoint_r_tmp = kpoint_r

    !This is to be sure that the k point fractional coordinates are stored.
    !If they are, the elec_read_band_energy will return.
!    call elec_read_band_energy

    ! We will call this only if we have not read in the cell before. With a parallel build
    ! the non root nodes would not get the updated cell and perform another bohr2ang
    ! conversion. If I try to insert a comms_bcast in the cell_get_
    ! F.Mildner 04/2023
    if (.not. abs(cell_volume) .gt. 0.0_dp) then
      call cell_get_real_lattice
      call cell_calc_lattice
    end if
    do loop = 1, num_kpoints_on_node(my_node_id)
      call utility_reciprocal_frac_to_cart(kpoint_r_tmp(:, loop), kpoint_r_cart_tmp(:, loop), recip_lattice)
!      print*,kpoint_r_tmp(1,loop),kpoint_r_tmp(2,loop),kpoint_r_tmp(3,loop),&
!       kpoint_r_cart_tmp(1,loop),kpoint_r_cart_tmp(2,loop),kpoint_r_cart_tmp(3,loop)
!print*,loop,kpoint_r_cart_tmp(:,loop),recip_lattice
    end do

    kpoint_r_cart = kpoint_r_cart_tmp
!      do loop=1,num_kpoints_on_node(my_node_id)
!      print*,kpoint_r_tmp(1,loop),kpoint_r_cart(1,loop),&
!      kpoint_r_tmp(2,loop),kpoint_r_cart(2,loop),&
!      kpoint_r_tmp(3,loop),kpoint_r_cart(3,loop)
!      end do

    deallocate (kpoint_r_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: cell_calc_kpoint_r_cart - &
&     failed to deallocate kpoint_r_tmp')

    deallocate (kpoint_r_cart_tmp, stat=ierr)
    if (ierr /= 0) call io_error('Error: cell_calc_kpoint_r_cart - &
&     failed to deallocate kpoint_r_cart_tmp')

  end subroutine cell_calc_kpoint_r_cart
end module od_cell
