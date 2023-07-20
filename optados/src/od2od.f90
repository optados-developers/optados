!-*- mode: F90 -*-!
module od_conv
  !! Helper module for od2od used for file conversions.

  use od_constants, only: dp
  use od_electronic, only: elec_read_optical_mat, elec_read_band_gradient, elec_read_elnes_mat,&
       & elec_pdos_read, elec_read_band_energy, omefile_header, domefile_header, pdosfile_header,&
       & elnesfile_header, elec_read_foptical_mat, femfile_header
  use od_parameters, only: iprint
  use od_io, only: stdout, io_error, seedname
  implicit none

  character(len=80), save :: outseedname
  !! It's conceivable that you might not
  !!  want to write over what you already have.
  character(len=80), save ::  infile
  !! Type of file to convert from.
  character(len=80), save :: outfile
  !! Type of file to convert to.
  character(len=10), save :: format_precision = "es23.10"
  !! Things get messy below 10 s.f. between bin files and fmt files
contains
  !=========================================================================
  subroutine print_usage()
    !! Writes the usage of the program to stdout
    write (stdout, '(A)')
    write (stdout, '(A)') " OptaDOS od2od "
    write (stdout, '(A)')
    write (stdout, '(A)') " Usage: od2od -i/--in_file <in_type> -o/--out_file <out_type> -w/--out_seedname [seedout] [seedname] "
    write (stdout, '(A)')
    write (stdout, '(A)') " [seedname] and [seedout] are optional input and output seednames"
    write (stdout, '(A)')
    write (stdout, '(A)') " <in_type> and <out_type> is one of: "
    write (stdout, '(A)') "       ome_fmt : a formatted optical matrix element file"
    write (stdout, '(A)') "       ome_bin : an unformatted optical matrix element file"
    ! Added by F. Mildner (04/2023) for photoemission
    write (stdout, '(A)') "       fem_fmt : a formatted free electron optical matrix element file"
    write (stdout, '(A)') "       fem_bin : an unformatted free electron optical matrix element file"

    write (stdout, '(A)') "      dome_fmt : a formatted diagonal optical matrix element file"
    write (stdout, '(A)') "      dome_bin : an unformatted diagonal optical matrix element file"
    write (stdout, '(A)') "      pdos_fmt : a formatted projected density of states file"
    write (stdout, '(A)') "      pdos_bin : an unformatted projected density of states file"
    write (stdout, '(A)') "     elnes_fmt : an formatted ELNES file"
    write (stdout, '(A)') "     elnes_bin : an unformatted ELNES file"
    write (stdout, '(A)') "         dummy : no input or output file (for testing)"
    write (stdout, '(A)')
    write (stdout, '(A)') " Known issues: (1) a seedname.bands file also needs to be present until"
    write (stdout, '(A)') "    I've thought of a better way to do it."
    write (stdout, '(A)') "      (3) It only works in serial."
    write (stdout, '(A)') "      (4) It only decides if the output format is sane after it's"
    write (stdout, '(A)') "     read the input."
    write (stdout, '(A)') "      (5) I need to think more about the amount of precision in in a formatted"
    write (stdout, '(A)') "     out file."
    write (stdout, '(A)') "      (6) File versions and headers could be better stored and reproduced."
    write (stdout, '(A)')
    write (stdout, '(A)') " Features: (1) Ability to convert a ome into a dome."
    write (stdout, '(A)')
  end subroutine print_usage

  !=========================================================================
  subroutine conv_get_seedname
    !! Set the seedname from the command line
    use od_io, only: seedname
    implicit none

    integer :: num_arg
    integer :: i !! Temporary variable
    character(len=50) :: ctemp

    num_arg = command_argument_count()

    outseedname = 'optados'
    seedname = 'optados' !! set to optados until proven otherwise
    i = 1
    do while (i .le. num_arg)
      call get_command_argument(i, ctemp)
      select case (trim(ctemp))
      case ("-i", "--in_file")
        i = i + 1
        call get_command_argument(i, infile)
      case ("-o", "--out_file")
        i = i + 1
        call get_command_argument(i, outfile)
      case ("-w", "--out_seedname")
        i = i + 1
        call get_command_argument(i, outseedname)
      case ("--") !! End of flags
        i = i + 1
        call get_command_argument(i, seedname)
      case default
        if (seedname == 'optados') then
          seedname = trim(ctemp)
          i = i + 1
        else !! We've already set the seedname so it can't be that again!
          call print_usage
          call io_error('Wrong command line arguments, see logfile for usage')
        end if
      end select
      i = i + 1
    end do

    if (outseedname == 'optados') then
      outseedname = seedname
    end if

  end subroutine conv_get_seedname

  !=========================================================================
  ! O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================

  !=========================================================================
  subroutine read_ome_fmt()
    !! Read a formatted Optical Matrix Elements file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, optical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string, string2
    integer :: ik, is, ib, i, jb, ome_unit = 6

    write (stdout, *) " Read a formatted ome file. "

    if (.not. allocated(optical_mat)) then
      write (stdout, *) " Allocating optical_mat."
      allocate (optical_mat(nbands, nbands, 3, nkpoints, nspins))
    end if

    open (unit=ome_unit, form='formatted', recl=1073741824, file=trim(seedname)//".ome_fmt")

    ! Total number of elements of ome
    write (string, '(I0,"(1x,",a,")")') 3*nbands*nbands, trim(format_precision)
    ! write(stdout,*) string

    ! write(string,'(a)') trim(format_precision)

    read (ome_unit, '('//trim(format_precision)//')') file_version

    read (ome_unit, '(a80)') omefile_header

    ! write(0,*) nkpoints, nspins, nbands

    do ik = 1, nkpoints
      do is = 1, nspins
        read (ome_unit, '('//trim(string)//')') (((optical_mat(ib, jb, i, ik, is), ib=1, nbands), &
             &jb=1, nbands), i=1, 3)
      end do
    end do

    optical_mat = optical_mat*(bohr2ang*H2eV)

    close (unit=ome_unit)

    write (stdout, *) trim(seedname)//".ome_fmt"//"--> Formatted ome sucessfully read. "

  end subroutine read_ome_fmt

  !=========================================================================
  subroutine write_ome_fmt()
    !! Write a formatted ome file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, optical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, ome_unit = 6

    write (stdout, *) " Write a formatted ome file. "

    optical_mat = optical_mat/(bohr2ang*H2eV)

    open (unit=ome_unit, form='formatted', file=trim(outseedname)//".ome_fmt")

    write (string, '(I0,"(1x,",a,")")') 3*nbands*nbands, trim(format_precision)
    !   write(stdout,*) string

    write (stdout, '(a80)') omefile_header
    write (stdout, '(a80)') adjustl(omefile_header)

    write (ome_unit, '('//trim(format_precision)//')') file_version
    write (ome_unit, '(a80)') adjustl(omefile_header)

    do ik = 1, nkpoints
      do is = 1, nspins
        write (ome_unit, '('//trim(string)//')') (((optical_mat(ib, jb, i, ik, is), ib=1, nbands), &
             &jb=1, nbands), i=1, 3)
      end do
    end do

    close (unit=ome_unit)

    write (stdout, *) " Sucesfully written a formatted ome file --> "//trim(outseedname)//".ome_fmt"
  end subroutine write_ome_fmt

  !=========================================================================
  subroutine read_ome_bin()
    !! Read a binary ome file. Wrapper to keep the naming tidy.
    implicit none
    write (stdout, *) " Read a formatted ome file. "

    call elec_read_optical_mat()
    write (stdout, *) " "//trim(seedname)//".ome_bin"//"--> Unformatted ome sucessfully read. "
  end subroutine read_ome_bin

  !=========================================================================
  subroutine write_ome_bin()
    !! Write a binary ome file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, optical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, ome_unit = 6

    write (stdout, *) " Write a binary ome file."

    optical_mat = optical_mat/(bohr2ang*H2eV)

    open (unit=ome_unit, form='unformatted', file=trim(outseedname)//".ome_bin")

    write (stdout, *) "-> Omefile_version ", file_version
    write (ome_unit) file_version
    write (stdout, *) "-> Omefile_header ", trim(omefile_header)
    write (ome_unit) adjustl(omefile_header)

    ! write(0,*) nkpoints, nspins, nbands
    do ik = 1, nkpoints
      do is = 1, nspins
        write (ome_unit) (((optical_mat(ib, jb, i, ik, is), ib=1, nbands), &
             &jb=1, nbands), i=1, 3)
      end do
    end do

    write (stdout, *) " Sucesfully written an unformatted ome file --> "//trim(outseedname)//".ome_bin"
  end subroutine write_ome_bin

  !=========================================================================
  ! F R E E   E L E C T R O N   O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================

  !=========================================================================
  subroutine read_fem_fmt()
    !! Read a formatted Optical Matrix Elements file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, foptical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string, string2
    integer :: ik, is, ib, i, jb, fem_unit = 6

    write (stdout, *) " Read a formatted .fem file. "

    if (.not. allocated(foptical_mat)) then
      write (stdout, *) " Allocating foptical_mat."
      allocate (foptical_mat(nbands + 1, nbands + 1, 3, nkpoints, nspins))
    end if

    open (unit=fem_unit, form='formatted', recl=1073741824, file=trim(seedname)//".fem_fmt")

    ! Total number of elements of ome
    write (string, '(I0,"(1x,",a,")")') 3*(nbands + 1)*(nbands + 1), trim(format_precision)
    ! write(stdout,*) string

    ! write(string,'(a)') trim(format_precision)

    read (fem_unit, '('//trim(format_precision)//')') file_version

    read (fem_unit, '(a80)') omefile_header

    ! write(0,*) nkpoints, nspins, nbands

    do ik = 1, nkpoints
      do is = 1, nspins
        read (fem_unit, '('//trim(string)//')') (((foptical_mat(ib, jb, i, ik, is), ib=1, nbands + 1), &
             &jb=1, nbands + 1), i=1, 3)
      end do
    end do

    foptical_mat = foptical_mat*(bohr2ang*H2eV)

    close (unit=fem_unit)

    write (stdout, *) trim(seedname)//".fem_fmt"//"--> Formatted fem sucessfully read. "

  end subroutine read_fem_fmt

  !=========================================================================
  subroutine write_fem_fmt()
    !! Write a formatted ome file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, foptical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, fem_unit = 6

    write (stdout, *) " Write a formatted .fem file. "

    foptical_mat = foptical_mat/(bohr2ang*H2eV)

    open (unit=fem_unit, form='formatted', file=trim(outseedname)//".fem_fmt")

    write (string, '(I0,"(1x,",a,")")') 3*(nbands + 1)*(nbands + 1), trim(format_precision)
    !   write(stdout,*) string

    write (stdout, '(a80)') femfile_header
    write (stdout, '(a80)') adjustl(femfile_header)

    write (fem_unit, '('//trim(format_precision)//')') file_version
    write (fem_unit, '(a80)') adjustl(femfile_header)

    do ik = 1, nkpoints
      do is = 1, nspins
        write (fem_unit, '('//trim(string)//')') (((foptical_mat(ib, jb, i, ik, is), ib=1, nbands + 1), &
             &jb=1, nbands + 1), i=1, 3)
      end do
    end do

    close (unit=fem_unit)

    write (stdout, *) " Sucesfully written a formatted fem file --> "//trim(outseedname)//".fem_fmt"
  end subroutine write_fem_fmt

  !=========================================================================
  subroutine read_fem_bin()
    !! Read a binary ome file. Wrapper to keep the naming tidy.
    implicit none
    write (stdout, *) " Read a formatted ome file. "

    call elec_read_foptical_mat()
    write (stdout, *) " "//trim(seedname)//".fem_bin"//"--> Unformatted ome sucessfully read. "
  end subroutine read_fem_bin

  !=========================================================================
  subroutine write_fem_bin()
    !! Write a binary ome file.
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, foptical_mat
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, fem_unit = 6

    write (stdout, *) " Write a binary fem file."

    foptical_mat = foptical_mat/(bohr2ang*H2eV)

    open (unit=fem_unit, form='unformatted', file=trim(outseedname)//".fem_bin")

    write (stdout, *) "-> Femfile_version ", file_version
    write (fem_unit) file_version
    write (stdout, *) "-> Femfile_header ", trim(femfile_header)
    write (fem_unit) adjustl(femfile_header)

    ! write(0,*) nkpoints, nspins, nbands
    do ik = 1, nkpoints
      do is = 1, nspins
        write (fem_unit) (((foptical_mat(ib, jb, i, ik, is), ib=1, nbands + 1), &
             &jb=1, nbands + 1), i=1, 3)
      end do
    end do

    write (stdout, *) " Sucesfully written an unformatted fem file --> "//trim(outseedname)//".fem_bin"
  end subroutine write_fem_bin

  !=========================================================================
  ! D I A G O N A L  O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================

  !=========================================================================
  subroutine read_dome_fmt()
    !! Read a diagonal ome formatted file.
    use od_constants, only: dp
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, band_gradient
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, dome_unit = 6

    write (stdout, *) " Read a formatted dome file. "

    if (.not. allocated(band_gradient)) then
      write (stdout, *) " Allocating band_gradient"
      allocate (band_gradient(nbands, 3, nkpoints, nspins))
    end if

    open (unit=dome_unit, form='formatted', file=trim(seedname)//".dome_fmt")

    write (string, '(i0,"(1x,",a,")")') 3*nbands, trim(format_precision)

    read (dome_unit, '('//trim(format_precision)//')') file_version

    read (dome_unit, '(a80)') domefile_header

    do ik = 1, nkpoints
      do is = 1, nspins
        read (dome_unit, '('//trim(string)//')') ((band_gradient(ib, i, ik, is), ib=1, nbands), &
             &i=1, 3)
      end do
    end do

    band_gradient = band_gradient*(bohr2ang*H2eV)

    close (unit=dome_unit)

    write (stdout, *) " "//trim(seedname)//".dome_fmt"//"--> Formatted ome sucessfully read. "
  end subroutine read_dome_fmt

  !=========================================================================
  subroutine write_dome_fmt()
    !! Write a diagonal ome formatted file.
    use od_constants, only: dp
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, band_gradient
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version
    character(len=100):: string
    integer :: ik, is, ib, i, jb, dome_unit = 6

    write (stdout, *) " Write a formatted ome file."

    open (unit=dome_unit, form='formatted', file=trim(outseedname)//".dome_fmt")

    write (string, '(I0,"(1x,",a,")")') 3*nbands, trim(format_precision)

    write (dome_unit, '('//trim(format_precision)//')') file_version
    write (dome_unit, '(a80)') adjustl(domefile_header)

    band_gradient = band_gradient/(bohr2ang*H2eV)

    do ik = 1, nkpoints
      do is = 1, nspins
        write (dome_unit, '('//trim(string)//')') ((band_gradient(ib, i, ik, is), ib=1, nbands), &
                                                   i=1, 3)
      end do
    end do

    close (unit=dome_unit)

    write (stdout, *) " Sucesfully written a formatted dome file --> "//trim(outseedname)//".dome_fmt"
  end subroutine write_dome_fmt

  !=========================================================================
  subroutine read_dome_bin()
    !! Read a diagonal ome file. Wrapper to keep the naming scheme tidy.
    implicit none
    write (stdout, *) " Read a binary dome file."
    call elec_read_band_gradient()
  end subroutine read_dome_bin

  !=========================================================================
  subroutine write_dome_bin()
    !! Write a diagonal ome file.
    use od_constants, only: dp
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, band_gradient
    use od_constants, only: bohr2ang, H2eV
    implicit none

    real(dp):: file_version = 1.0_dp          ! File version

    integer :: ik, is, ib, i, jb, dome_unit = 6

    write (stdout, *) " Write a binary dome file."

    open (unit=dome_unit, form='unformatted', file=trim(outseedname)//".dome_bin")

    band_gradient = band_gradient/(bohr2ang*H2eV)

    write (dome_unit) file_version
    write (dome_unit) adjustl(domefile_header)

    do ik = 1, nkpoints
      do is = 1, nspins
        write (dome_unit) ((band_gradient(ib, i, ik, is), ib=1, nbands), &
             &i=1, 3)
      end do
    end do
    write (stdout, *) " Sucesfully written a binary dome file --> "//trim(outseedname)//".dome_bin"
  end subroutine write_dome_bin

  !=========================================================================
  ! P R O J E C T E D   D O S
  !=========================================================================

  !=========================================================================
  subroutine write_pdos_fmt()
    !! Write a formatted pdos file.
    use od_electronic, only: pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights
    use od_cell, only: nkpoints, kpoint_r
    use od_io, only: stdout, seedname, io_file_unit
    implicit none

    integer :: ik, is, ib
    integer :: pdos_in_unit
    character(len=80) :: string, string2
    real(dp) :: file_version = 1.0_dp

    write (stdout, *) " Write a formatted pdos file."
    !-------------------------------------------------------------------------!
    ! W R I T E   T H E   D A T A   H E A D E R

    pdos_in_unit = io_file_unit()

    open (unit=pdos_in_unit, file=trim(outseedname)//".pdos_fmt", form='formatted')
    write (pdos_in_unit, '('//trim(format_precision)//')') file_version
    write (pdos_in_unit, '(a80)') adjustl(pdosfile_header)

    write (pdos_in_unit, '(a10, i6)') "Kpoints", pdos_mwab%nkpoints
    write (pdos_in_unit, '(a10, i6)') "Spins", pdos_mwab%nspins
    write (pdos_in_unit, '(a10, i6)') "Orbials", pdos_mwab%norbitals
    write (pdos_in_unit, '(a10, i6)') "Bands", pdos_mwab%nbands

    !write(stdout,'(a30,i6)') "DEBUG: pdos_mwab%nkpoints= ",pdos_mwab%nkpoints
    !write(stdout,'(a30,i6)') "DEBUG: pdos_mwab%nspins= ",pdos_mwab%nspins
    !write(stdout,'(a30,i6)') "DEBUG: pdos_mwab%norbitals= ",pdos_mwab%norbitals
    !write(stdout,'(a30,i6)') "DEBUG: pdos_mwab%nbands= ",pdos_mwab%nbands

    ! These should all be allocated!
    !allocate(pdos_orbital%species_no(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    !allocate(pdos_orbital%rank_in_species(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    !allocate(pdos_orbital%am_channel(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")

    write (string, '(i7,"(1x,",a,")")') pdos_mwab%norbitals, "i5"

    write (string2, '(i7,"(1x,",a,")")') pdos_mwab%norbitals, trim(format_precision)

    write (pdos_in_unit, '(a60)') " Species number for each orbital"
    write (pdos_in_unit, '('//trim(string)//')') pdos_orbital%species_no(1:pdos_mwab%norbitals)
    write (pdos_in_unit, '(a60)') " Species rank for each orbital"
    write (pdos_in_unit, '('//trim(string)//')') pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    write (pdos_in_unit, '(a60)') " AM channel for each orbital"
    write (pdos_in_unit, '('//trim(string)//')') pdos_orbital%am_channel(1:pdos_mwab%norbitals)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A

    ! These should already be alloacted
    ! allocate(nbands_occ(1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate nbands_occ"
    ! allocate(pdos_weights(1:pdos_mwab%norbitals,1:pdos_mwab%nbands, &
    !     1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate pdos_weights"

    do ik = 1, nkpoints
      ! The kpoint number, followed by the kpoint-vector
      write (pdos_in_unit, '(i6,3'//trim(format_precision)//')') ik, kpoint_r(:, ik)
      do is = 1, pdos_mwab%nspins
        write (pdos_in_unit, '(i6)') is ! this is the spin number
        write (pdos_in_unit, '(i6)') nbands_occ(ik, is)
        do ib = 1, nbands_occ(ik, is)

          !      write(stdout,*) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
          !      write(stdout,*) " DEBUG:", ib, ik, is
          !      write(stdout,*) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "

          write (pdos_in_unit, '('//trim(string2)//')') pdos_weights(1:pdos_mwab%norbitals, ib, ik, is)
        end do
      end do
    end do

    close (pdos_in_unit)

    write (stdout, *) " Sucesfully written a formtted pdos file --> "//trim(outseedname)//".dome_bin"

  end subroutine write_pdos_fmt

  !=========================================================================
  subroutine read_pdos_fmt()
    !! Read a formatted pdos file.
    use od_electronic, only: pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights
    use od_cell, only: nkpoints, kpoint_r, num_kpoints_on_node
    use od_io, only: stdout, seedname, io_file_unit
    use od_comms, only: my_node_id
    implicit none

    integer :: ik, is, ib, idummy, ierr, io
    integer :: pdos_in_unit
    character(len=80) :: string, dummy, string2
    real(dp) :: file_version

    write (stdout, *) " Read a formatted pdos file."
    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R

    file_version = 1.0_dp
    pdos_in_unit = io_file_unit()

    open (unit=pdos_in_unit, file=trim(seedname)//".pdos_fmt", form='formatted')
    read (pdos_in_unit, '('//trim(format_precision)//')') file_version
    read (pdos_in_unit, '(a80)') pdosfile_header

    read (pdos_in_unit, '(a10, i6)') dummy, pdos_mwab%nkpoints
    read (pdos_in_unit, '(a10, i6)') dummy, pdos_mwab%nspins
    read (pdos_in_unit, '(a10, i6)') dummy, pdos_mwab%norbitals
    read (pdos_in_unit, '(a10, i6)') dummy, pdos_mwab%nbands

    !write(stdout,'(a, i6)') "DEBUG: pdos_mwab%nkpoints= ",pdos_mwab%nkpoints
    !write(stdout,'(a, i6)') "DEBUG: pdos_mwab%nspins= ",pdos_mwab%nspins
    !write(stdout,'(a, i6)') "DEBUG: pdos_mwab%norbitals= ",pdos_mwab%norbitals
    !write(stdout,'(a, i6)') "DEBUG: pdos_mwab%nbands= ",pdos_mwab%nbands

    write (string, '(i7,"(1x,",a,")")') pdos_mwab%norbitals, "i5"

    write (string2, '(i7,"(1x,",a,")")') pdos_mwab%norbitals, trim(format_precision)

    ! These should all be allocated!
    allocate (pdos_orbital%species_no(pdos_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
    allocate (pdos_orbital%rank_in_species(pdos_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")
    allocate (pdos_orbital%am_channel(pdos_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(" Error : cannot allocate pdos_orbital")

    allocate (nbands_occ(1:num_kpoints_on_node(my_node_id), 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate nbands_occ"

    allocate (pdos_weights(1:pdos_mwab%norbitals, 1:pdos_mwab%nbands, &
                           1:num_kpoints_on_node(my_node_id), 1:pdos_mwab%nspins), stat=ierr)
    if (ierr /= 0) stop " Error : cannot allocate pdos_weights"

    read (pdos_in_unit, '(a60)') dummy
    read (pdos_in_unit, '('//trim(string)//')') pdos_orbital%species_no(1:pdos_mwab%norbitals)
    read (pdos_in_unit, '(a60)') dummy
    read (pdos_in_unit, '('//trim(string)//')') pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    read (pdos_in_unit, '(a60)') dummy
    read (pdos_in_unit, '('//trim(string)//')') pdos_orbital%am_channel(1:pdos_mwab%norbitals)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A

    ! These should already be alloacted
    ! allocate(nbands_occ(1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate nbands_occ"
    ! allocate(pdos_weights(1:pdos_mwab%norbitals,1:pdos_mwab%nbands, &
    !     1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate pdos_weights"

    do ik = 1, nkpoints
      ! The kpoint number, followed by the kpoint-vector
      read (pdos_in_unit, '(i6,3'//trim(format_precision)//')') idummy, kpoint_r(:, ik)
      do is = 1, pdos_mwab%nspins
        read (pdos_in_unit, '(i6)') idummy ! this is the spin number
        read (pdos_in_unit, '(i6)') nbands_occ(ik, is)
        do ib = 1, nbands_occ(ik, is)

          !      write(stdout,*) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
          !      write(stdout,*) " DEBUG:", ib, ik, is
          !      write(stdout,*) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "

          read (pdos_in_unit, '('//trim(string2)//')') (pdos_weights(io, ib, ik, is), io=1, pdos_mwab%norbitals)
        end do
      end do
    end do

    close (pdos_in_unit)

    write (stdout, *) " "//trim(seedname)//".pdos_fmt"//"--> Formatted pdos sucessfully read. "

  end subroutine read_pdos_fmt

  !=========================================================================
  subroutine read_pdos_bin()
    !! Wrapper to read a pdos binary file. Useful to keep the code tidy.
    implicit none
    write (stdout, *) " Read a binary pdos file."
    call elec_pdos_read()
  end subroutine read_pdos_bin

  !=========================================================================
  subroutine write_pdos_bin()
    !! Write a binary pdos file
    use od_constants, only: dp
    use od_electronic, only: pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights
    use od_cell, only: nkpoints, kpoint_r
    use od_io, only: io_file_unit

    implicit none

    !  integer,parameter:: num_popn_orb=1          ! Number of LCAO projectors
    ! integer,parameter:: max_eigenv=1       ! Number of bands included in matrix elements
    real(dp):: file_version = 1.0_dp  ! File version
    ! integer:: species(1:num_popn_orb)! Atomic species associated with each projector
    ! integer:: ion(1:num_popn_orb)    ! Ion associated with each projector
    ! integer:: am_channel(1:num_popn_orb)     ! Angular momentum channel
    ! integer:: num_eigenvalues(1:nspins)   ! Number of eigenvalues per spin channel
    ! real(dp):: kpoint_positions(1:nkpoints,1:3) ! k_x, k_y, k_z in fractions of BZ
    ! real(dp):: pdos_weights(1:num_popn_orb,max_eigenv,num_kpoints,num_spins)!Matrix elements
    !character(len=80):: file_header ! File header comment

    integer :: ik, is, ib, pdos_file, io, idex

    pdos_file = io_file_unit()

    open (unit=pdos_file, file=trim(outseedname)//".pdos_bin", form='unformatted')

    write (stdout, *) " Write a binary pdos file."

    write (pdos_file) file_version
    write (pdos_file) adjustl(pdosfile_header)
    write (pdos_file) pdos_mwab%nkpoints
    write (pdos_file) pdos_mwab%nspins
    write (pdos_file) pdos_mwab%norbitals
    write (pdos_file) pdos_mwab%nbands
    write (pdos_file) pdos_orbital%species_no(1:pdos_mwab%norbitals)
    write (pdos_file) pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    write (pdos_file) pdos_orbital%am_channel(1:pdos_mwab%norbitals)

    write (stdout, *) pdos_mwab%nkpoints
    write (stdout, *) pdos_mwab%nspins
    write (stdout, *) pdos_mwab%norbitals
    write (stdout, *) pdos_mwab%nbands
    write (stdout, *) pdos_orbital%species_no(1:pdos_mwab%norbitals)
    write (stdout, *) pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    write (stdout, *) pdos_orbital%am_channel(1:pdos_mwab%norbitals)

    do ik = 1, pdos_mwab%nkpoints
      write (stdout, *) "loop", ik
      write (pdos_file) ik, (kpoint_r(idex, ik), idex=1, 3)
      do is = 1, pdos_mwab%nspins
        write (pdos_file) is
        write (pdos_file) nbands_occ(ik, is)
        write (stdout, *) is, nbands_occ(ik, is)
        do ib = 1, nbands_occ(ik, is)
          write (pdos_file) (pdos_weights(io, ib, ik, is), io=1, pdos_mwab%norbitals)
        end do
      end do
    end do

    write (stdout, *) " Sucesfully written a binary pdos file --> "//trim(outseedname)//".pdos_bin"

  end subroutine write_pdos_bin

  !=========================================================================
  ! E L N E S   M A T R I X   E L E M E N T S
  !=========================================================================

  !=========================================================================
  subroutine read_elnes_fmt()
    !! Read a formatted elnes file.
    use od_electronic, only: elec_elnes_find_channel_names, elnes_orbital, &
         & elnes_mwab, elnes_mat
    use od_io, only: io_file_unit, seedname
    use od_cell, only: num_kpoints_on_node
    implicit none

    character(len=20) :: dummy20, dummy10

    real(dp) :: file_version
    !! The file verioson format to write. Currently we're on version 1.
    integer :: elnes_unit, ik, is, iorb, ib, indx, ierr
    !! Loop variable.
    character(len=80) :: string, string2
    !! Tempory string manipulation variable.

    write (stdout, *) " Read a formatted elnes file."

    elnes_unit = io_file_unit()

    open (unit=elnes_unit, file=trim(seedname)//".elnes_fmt", form='formatted')

    read (elnes_unit, '('//trim(format_precision)//')') file_version
    read (elnes_unit, '(a80)') elnesfile_header

    read (elnes_unit, '(a20,1x,i5)') dummy20, elnes_mwab%norbitals
    read (elnes_unit, '(a20,1x,i5)') dummy20, elnes_mwab%nbands
    read (elnes_unit, '(a20,1x,i5)') dummy20, elnes_mwab%nkpoints
    read (elnes_unit, '(a20,1x,i5)') dummy20, elnes_mwab%nspins

    write (string, '(i7,"(1x,",a,")")') elnes_mwab%norbitals, "i5"
    write (string2, '(i7,"(1x,",a,")")') elnes_mwab%norbitals*elnes_mwab%nbands*3, trim(format_precision)

    allocate (elnes_orbital%ion_no(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbital%ion_no')
    allocate (elnes_orbital%species_no(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbital%species_no')
    allocate (elnes_orbital%rank_in_species(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbitall%rank_in_species')
    allocate (elnes_orbital%shell(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbitall%shell')
    allocate (elnes_orbital%am_channel(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbital%am_channel')
    allocate (elnes_orbital%am_channel_name(elnes_mwab%norbitals), stat=ierr)
    if (ierr /= 0) call io_error(' Error : read_elnes_fmt cannot allocate elnes_orbital%am_channel_name')

    read (elnes_unit, '(a10,'//trim(string)//')') dummy10, elnes_orbital%species_no(1:elnes_mwab%norbitals)
    read (elnes_unit, '(a10,'//trim(string)//')') dummy10, elnes_orbital%rank_in_species(1:elnes_mwab%norbitals)
    read (elnes_unit, '(a10,'//trim(string)//')') dummy10, elnes_orbital%shell(1:elnes_mwab%norbitals)
    read (elnes_unit, '(a10,'//trim(string)//')') dummy10, elnes_orbital%am_channel(1:elnes_mwab%norbitals)

    allocate (elnes_mat(1:elnes_mwab%norbitals, 1:elnes_mwab%nbands, 1:3, &
                        1:num_kpoints_on_node(0), 1:elnes_mwab%nspins), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating elnes_mat in read_elnes_fmt')

    do ik = 1, num_kpoints_on_node(0)
      do is = 1, elnes_mwab%nspins
        read (elnes_unit, '('//trim(string2)//')') (((elnes_mat(iorb, ib, indx, ik, is), iorb=1, elnes_mwab%norbitals), &
                                                     ib=1, elnes_mwab%nbands), indx=1, 3)
      end do
    end do

    call elec_elnes_find_channel_names()

    write (stdout, *) " "//trim(seedname)//".elnes_fmt"//"--> Formatted elnes sucessfully read. "

  end subroutine read_elnes_fmt

  !=========================================================================
  subroutine write_elnes_fmt()
    !! Soubroute to write a formatted elnes file.
    use od_electronic, only: elnes_mwab, elnes_orbital, elnes_mat, &
         & elec_elnes_find_channel_numbers
    use od_cell, only: num_kpoints_on_node
    use od_io, only: io_file_unit
    implicit none

    real(dp) :: file_version = 1.0_dp
    !! The file verioson format to write. Currently we're on version 1.
    integer :: elnes_unit
    !! File unit number to write to.
    integer :: ik, is, iorb, ib, indx
    !! Loop variables.
    character(len=80) :: string, string2
    !! Tempory string manipulation variables.

    write (stdout, *) " Write a formatted elnes file."

    ! CASTEP (hence the bin file) and OptaDOS think about am_channel numbers
    ! differently. To keep consistent we convert to CASTEP's numbering scheme
    ! before we write out.
    call elec_elnes_find_channel_numbers()

    elnes_unit = io_file_unit()

    open (unit=elnes_unit, file=trim(outseedname)//".elnes_fmt", form='formatted')

    write (elnes_unit, '('//trim(format_precision)//')') file_version
    write (elnes_unit, '(a80)') adjustl(elnesfile_header)

    write (elnes_unit, '(a20,1x,i5)') "Norbitals", elnes_mwab%norbitals
    write (elnes_unit, '(a20,1x,i5)') "Nbands", elnes_mwab%nbands
    write (elnes_unit, '(a20,1x,i5)') "Nkpoints", elnes_mwab%nkpoints
    write (elnes_unit, '(a20,1x,i5)') "Nspins", elnes_mwab%nspins

    write (string, '(i7,"(1x,",a,")")') elnes_mwab%norbitals, "i5"
    write (string2, '(i7,"(1x,",a,")")') elnes_mwab%norbitals*elnes_mwab%nbands*3, trim(format_precision)

    write (elnes_unit, '(a10,'//trim(string)//')') "Species_no", elnes_orbital%species_no(1:elnes_mwab%norbitals)
    write (elnes_unit, '(a10,'//trim(string)//')') "Rank", elnes_orbital%rank_in_species(1:elnes_mwab%norbitals)
    write (elnes_unit, '(a10,'//trim(string)//')') "Shell", elnes_orbital%shell(1:elnes_mwab%norbitals)
    write (elnes_unit, '(a10,'//trim(string)//')') "Am_channel", elnes_orbital%am_channel(1:elnes_mwab%norbitals)

    do ik = 1, num_kpoints_on_node(0)
      do is = 1, elnes_mwab%nspins
        write (elnes_unit, '('//trim(string2)//')') (((elnes_mat(iorb, ib, indx, ik, is), iorb=1, elnes_mwab%norbitals), &
                                                      ib=1, elnes_mwab%nbands), indx=1, 3)
      end do
    end do

    close (elnes_unit)

    write (stdout, *) " Sucesfully written a formatted elnes file --> "//trim(outseedname)//".elnes_fmt"

  end subroutine write_elnes_fmt

  !=========================================================================
  subroutine read_elnes_bin()
    !! Wrapper to read a binary elnes file. The wrapping allows us to have
    !! consistent names within the module, which makes life easier.
    implicit none
    write (stdout, *) " Read a binary elnes file."
    call elec_read_elnes_mat()
  end subroutine read_elnes_bin

  !=========================================================================
  subroutine write_elnes_bin()
    !! Writes a binary elnes file.
    use od_electronic, only: elec_elnes_find_channel_numbers, elnes_orbital,&
         & elnes_mat, elnes_mwab
    use od_cell, only: num_kpoints_on_node
    use od_io, only: io_file_unit
    implicit none

    real(dp) :: file_version = 1.0_dp
    !! The file verioson format to write. Currently we're on version 1.
    integer :: ik, is, ib, iorb, indx
    !! Loop variables
    integer :: elnes_unit
    !! File unit number to write to.

    write (stdout, *) " Write a binary elnes file."

    !write !! Some headers here?

    ! CASTEP (hence the bin file) and OptaDOS think about am_channel numbers
    ! differently. To keep consistent we convert to CASTEP's numbering scheme
    ! before we write out.
    call elec_elnes_find_channel_numbers()

    elnes_unit = io_file_unit()

    open (unit=elnes_unit, file=trim(outseedname)//".elnes_bin", form='unformatted')

    write (elnes_unit) file_version
    write (elnes_unit) adjustl(elnesfile_header)

    write (elnes_unit) elnes_mwab%norbitals
    write (elnes_unit) elnes_mwab%nbands
    write (elnes_unit) elnes_mwab%nkpoints
    write (elnes_unit) elnes_mwab%nspins

    ! write(string,'(i7,"(1x,",a,")")') elnes_mwab%norbitals,"i5"
    ! write(string2,'(i7,"(1x,",a,")")') elnes_mwab%norbitals*elnes_mwab%nbands*3, trim(format_precision)

    write (elnes_unit) elnes_orbital%species_no(1:elnes_mwab%norbitals)
    write (elnes_unit) elnes_orbital%rank_in_species(1:elnes_mwab%norbitals)
    write (elnes_unit) elnes_orbital%shell(1:elnes_mwab%norbitals)
    write (elnes_unit) elnes_orbital%am_channel(1:elnes_mwab%norbitals)

    do ik = 1, num_kpoints_on_node(0)
      do is = 1, elnes_mwab%nspins
        write (elnes_unit) (((elnes_mat(iorb, ib, indx, ik, is), iorb=1, elnes_mwab%norbitals), &
                             ib=1, elnes_mwab%nbands), indx=1, 3)
      end do
    end do

    write (stdout, *) " Sucesfully written a binary elnes file --> "//trim(outseedname)//".elnes_bin"

    close (elnes_unit)

  end subroutine write_elnes_bin

  !=========================================================================
  subroutine slice_an_ome()
    !! This routine takes OptaDOS's internal representation of an ome and
    !! puts its diagonal into its internal representation of a dome.
    !! It's useful for testing.
    use od_constants, only: dp
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, band_gradient, optical_mat
    implicit none

    integer :: loop
    !! Loop variable

    write (stdout, *) " Slicing an ome into a dome."

    if (.not. allocated(band_gradient)) then
      write (stdout, *) " Allocating band_gradient"
      allocate (band_gradient(nbands, 3, nkpoints, nspins))
    end if

    do loop = 1, nbands
      band_gradient(loop, :, :, :) = real(optical_mat(loop, loop, :, :, :), dp)
    end do

  end subroutine slice_an_ome

  !=========================================================================
  subroutine pad_an_ome()
    !! This routine takes OptaDOS's internal representation of an dome file, and
    !! uses it to construct an ome file, where the off diagonal elements are
    !! padded with zeros.
    !! It might be helpful if we're reading in other codes' input files.
    !! It might also be useful for testing OptaDOS itself.
    use od_constants, only: dp
    use od_io, only: io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error
    use od_cell, only: num_kpoints_on_node, nkpoints
    use od_electronic, only: nspins, nbands, band_gradient, optical_mat
    implicit none

    integer :: loop
    !! Loop variable.

    write (stdout, *) " Padding a dome into a ome."

    if (.not. allocated(optical_mat)) then
      write (stdout, *) " Allocating optical_mat"
      allocate (optical_mat(nbands, nbands, 3, nkpoints, nspins))
    end if

    ! Pad with zeros.
    optical_mat = 0.0_dp

    do loop = 1, nbands
      optical_mat(loop, loop, :, :, :) = band_gradient(loop, :, :, :)
    end do

  end subroutine pad_an_ome

  !=========================================================================
  subroutine report_arraysize()
    !! Write to stdout some info on the size of the arrays we're using. These
    !! are normally found in the .bands file.
    use od_electronic, only: nspins, nbands
    use od_cell, only: nkpoints
    use od_io, only: stdout
    implicit none

    write (stdout, '(a40,i5)') " Number of kpoints : ", nkpoints
    write (stdout, '(a40,i5)') " Number of bands : ", nbands
    write (stdout, '(a40,i5)') " Number of spins : ", nspins

  end subroutine report_arraysize

  !=========================================================================
  subroutine get_band_energy()
    !! Read the band file info which is prerequisite to know about k-points
    !! bands etc.
    !! It would be nice for od2od to work this out for itself. But at least this
    !! way it is consistent.  The problem is that to convert, say a ome to a ome
    !! one also requires a .bands file.
    implicit none

    write (stdout, *)
    write (stdout, *) "+----------------------------- K-point information --------------------------+"

    call elec_read_band_energy()
    call report_arraysize()
    write (stdout, *) "+----------------------------------------------------------------------------+"

  end subroutine get_band_energy

  !=========================================================================
  subroutine write_read_file()
    !! Noddy routine to prettify output
    implicit none
    write (stdout, *)
    write (stdout, *) "+-------------------------------- Read File ---------------------------------+"

  end subroutine write_read_file

end module od_conv

program od2od
  !! Program to convert checkpoint files from formatted to unformmated
  !! and vice versa - useful for switching between computers.
  !!
  !! Plan is to use to convert outputs of other DFT programs to ones that
  !! OptaDOS can read.
  !!
  !! AJM 2019
  use od_constants, only: dp
  use od_io, only: io_file_unit, stdout, stderr, io_error, seedname, io_time, &
    & io_date
  use od_conv
  use od_comms, only: num_nodes, comms_setup, comms_end
  implicit none

  logical :: file_found
  logical :: ome_conv, fem_conv, dome_conv, pdos_conv, elnes_conv, dummy_conv
  !! Flags to stop people trying to, say, read in a pdos and write out an
  !! elnes. That's not going to end well.
  real(kind=dp) :: time0, time1
  !! Varaibles for measuring exectuion time.
  character(len=9) :: pos
  !! Status and position of .odo file
  character(len=9) :: stat
  !! Position of .odo file
  character(len=9) :: ctime
  !! Temp. time string
  character(len=11):: cdate
  !! Temp. date string

  time0 = io_time()

  iprint = 4

  call comms_setup

  call conv_get_seedname

  stderr = io_file_unit()
  open (unit=stderr, file=trim(seedname)//'.opt_err')
  call io_date(cdate, ctime)
  write (stderr, *) 'od2od: Execution started on ', cdate, ' at ', ctime

  stdout = io_file_unit()
  open (unit=stdout, file=trim(seedname)//'.log')
  !-------------------------------------------------------------------------!
  write (stdout, *)
  write (stdout, *) 'od2od: Execution started on ', cdate, ' at ', ctime
  write (stdout, *)
  write (stdout, *) "+============================================================================+ "
  write (stdout, *) "|                                                                            | "
  write (stdout, *) "|                         OO   DDD   222    OO   DDD                         | "
  write (stdout, *) "|                        O  O  D  D     2  O  O  D  D                        | "
  write (stdout, *) "|                        O  O  D  D   22   O  O  D  D                        | "
  write (stdout, *) "|                        O  O  D  D  2     O  O  D  D                        | "
  write (stdout, *) "|                         OO   DDD   2222   OO   DDD                         | "
  write (stdout, *) "|                                                                            | "
  write (stdout, *) "|                   For doing the odd thing to OptaDOS files                 | "
  write (stdout, *) "|                                                                            | "
  write (stdout, *) "|                      OptaDOS Developers Group 2019 (C)                     | "
  write (stdout, *) "|                             (But blame Andrew)                             | "
  write (stdout, *) "|                                                                            | "
  write (stdout, *) "+============================================================================+ "

  if (num_nodes /= 1) then
    call io_error('od2od can only be used in serial...')
  end if

  write (stdout, *)
  write (stdout, *) "+--------------------------------- JOB CONTROL ------------------------------+"
  write (stdout, '(a40,i5)') "Number of nodes : ", num_nodes
  write (stdout, '(a40,a)') "Convert from : ", trim(infile)
  write (stdout, '(a40,a)') "Convert to : ", trim(outfile)
  write (stdout, '(a40,a)') "Seedname : ", trim(seedname)
  write (stdout, '(a40,a)') "Output Seedname : ", trim(outseedname)
  write (stdout, *) "+----------------------------------------------------------------------------+"

  ome_conv = .false.
  fem_conv = .false.
  dome_conv = .false.
  pdos_conv = .false.
  elnes_conv = .false.
  dummy_conv = .false.

  ! Main case to decide what file format to read in.
  read_input:select case(trim(infile))
case ("ome_fmt")
  ome_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_ome_fmt()
case ("ome_bin")
  ome_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_ome_bin()
case ("fem_fmt")
  fem_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_fem_fmt()
case ("fem_bin")
  fem_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_fem_bin()
case ("dome_fmt")
  dome_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_dome_fmt()
case ("dome_bin")
  dome_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_dome_bin()
case ("pdos_fmt")
  pdos_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_pdos_fmt()
case ("pdos_bin")
  pdos_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_pdos_bin()
case ("elnes_fmt")
  elnes_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_elnes_fmt()
case ("elnes_bin")
  elnes_conv = .true.
  call get_band_energy()
  call write_read_file()
  call read_elnes_bin()
case ("dummy")
  dummy_conv = .true.
  call get_band_energy()
  call write_read_file()
  write (stdout, *) " Not reading any input file."
case default
  call io_error('Unknown Input File format speccified')
  end select read_input
  write (stdout, *) "+----------------------------------------------------------------------------+"
  write (stdout, *)
  write (stdout, *) "+------------------------------- Write File ---------------------------------+"
  ! Main case to decide what file format to write.
  write_output:select case(trim(outfile))
case ("ome_fmt")
  if (.not. (dome_conv .or. ome_conv)) call io_error(' Input format '//trim(infile)//' not compatible with output format '&
       &//trim(outfile))
  if (dome_conv) call pad_an_ome()
  call write_ome_fmt()
case ("ome_bin")
  if (.not. (dome_conv .or. ome_conv)) call io_error(' Input format '//trim(infile)//' not compatible with output format '&
       &//trim(outfile))
  if (dome_conv) call pad_an_ome()
  call write_ome_bin()
case ("fem_fmt")
  if (.not. (fem_conv)) call io_error(' Input format '//trim(infile)//' not compatible with output format '&
       &//trim(outfile))
  call write_fem_fmt()
case ("fem_bin")
  if (.not. (fem_conv)) call io_error(' Input format '//trim(infile)//' not compatible with output format '&
       &//trim(outfile))
  call write_fem_bin()
case ("dome_fmt")
  if (.not. (dome_conv .or. ome_conv)) call io_error(' Input format '//trim(infile)//&
       &' not compatible with output format '//trim(outfile))
  if (ome_conv) call slice_an_ome()
  call write_dome_fmt()
case ("dome_bin")
  if (.not. (dome_conv .or. ome_conv)) call io_error(' Input format '//trim(infile)//&
       &' not compatible with output format '//trim(outfile))
  if (ome_conv) call slice_an_ome()
  call write_dome_bin()
case ("pdos_fmt")
  if (.not. pdos_conv) call io_error(' Input format '//trim(infile)//' not compatible with output format '//trim(outfile))
  call write_pdos_fmt()
case ("pdos_bin")
  if (.not. pdos_conv) call io_error(' Input format '//trim(infile)//' not compatible with output format '//trim(outfile))
  call write_pdos_bin()
case ("elnes_fmt")
  if (.not. elnes_conv) call io_error(' Input format '//trim(infile)//' not compatible with output format '//trim(outfile))
  call write_elnes_fmt()
case ("elnes_bin")
  if (.not. elnes_conv) call io_error(' Input format '//trim(infile)//' not compatible with output format '//trim(outfile))
  call write_elnes_bin()
case ("dummy")
  write (stdout, *) " Not writing any output file."
  if (dummy_conv) then
    write (stdout, *)
    write (stdout, *) "                 Dummy in + dummy out  -- who's the dummy now ?"
  else
    write (stdout, *)
    write (stdout, *) "                No point in taking up disk space unnecessarily, eh ?"
  end if
case default
  call io_error('Unknown Output File format speccified')
  end select write_output

  call io_date(cdate, ctime)
  write (stdout, *) "+----------------------------------------------------------------------------+"

  time1 = io_time()

  write (stdout, '(1x,a40,f11.3,a)') 'Total runtime :', time1 - time0, ' (sec)'
  write (stdout, *)
  write (stdout, *) 'od2od: Execution complete on ', cdate, ' at ', ctime
  write (stdout, *)

  close (stdout)
  close (stderr, status='delete')

  call comms_end

!  close(unit=stdout,status='delete')

end program od2od
