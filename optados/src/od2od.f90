!-*- mode: F90 -*-!
module  od_conv
  use od_constants, only: dp
  use od_electronic, only: elec_read_optical_mat, elec_read_band_gradient,  elec_read_elnes_mat,&
       & elec_pdos_read,  elec_read_band_energy, omefile_header, domefile_header, pdosfile_header,&
       & elnesfile_header
  use od_parameters, only: iprint
  use od_io, only:  stdout, io_error, seedname
  implicit none

  character(len=80), save :: outseedname ! It's conceivable that you might not
  !  want to write over what you already have.
  character(len=80), save ::  infile, outfile
 ! Things get messy below 10 s.f. between bin files and fmt files
  character(len=10), save :: format_precision="es23.10"
 
 contains
 !=========================================================================
   subroutine print_usage()
 !=========================================================================
     !! Writes the usage of the program to stdout
     write (stdout, '(A)')
     write (stdout, '(A)') " OptaDOS od2od "
     write (stdout, '(A)')
     write (stdout, '(A)') " Usage: od2od <in_type> <out_type> [seedname] [seedout]"
     write (stdout, '(A)')
     write (stdout, '(A)') " [seedname] and [seedout] are optional input and output seednames"
     write (stdout, '(A)')
     write (stdout, '(A)') " <in_type> and <out_type> is one of: "
     write (stdout, '(A)') "       ome_fmt : a formatted optical matrix element file"
     write (stdout, '(A)') "       ome_bin : an unformatted optical matrix element file"
     write (stdout, '(A)') "      dome_fmt : a formatted diagonal optical matrix element file"
     write (stdout, '(A)') "      dome_bin : an unformatted diagonal optical matrix element file"
     write (stdout, '(A)') "      pdos_fmt : a formatted projected density of states file"
     write (stdout, '(A)') "      pdos_bin : an unformatted projected density of states file"
     write (stdout, '(A)') "     elnes_fmt : an formatted ELNES file"
     write (stdout, '(A)') "     elnse_bin : an unformatted ELNES file"
     write (stdout, '(A)')
     write (stdout, '(A)') " Known issues: (1) a seedname.bands file also needs to be present until"
     write (stdout, '(A)') "    I've thought of a better way to do it."
     write (stdout, '(A)') "               (3) It only works in serial."
     write (stdout, '(A)') "               (4) It only decides if the output format is sane"
     write (stdout, '(A)') "    after it's read the input."
     write (stdout, '(A)') "               (5) I need to think more about the amount of precision in"
     write (stdout, '(A)') "    in a formatted out file."
     write (stdout, '(A)') "               (6) File versions and headers could be better stored and"
     write (stdout, '(A)') "    reproduced."
     write (stdout, '(A)')
     write (stdout, '(A)') " Features: (1) Ability to convert a ome into a dome."
     write (stdout, '(A)') 
  end subroutine print_usage 

  !=========================================================================
  subroutine conv_get_seedname
    !=========================================================================
    !! Set the seedname from the command line
    implicit none
    
    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 2) then
       seedname = 'optados'
    elseif (num_arg == 3) then
       call get_command_argument(3, seedname)
       outseedname=trim(seedname)
    elseif (num_arg == 4) then
       call get_command_argument(3, seedname)
       call get_command_argument(4, outseedname)
    else
       call print_usage
       call io_error('Wrong command line arguments, see logfile for usage')
    end if
    
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
   
 end subroutine conv_get_seedname

  !=========================================================================
  ! O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================

  
  !=========================================================================
  subroutine read_ome_fmt()
    !=========================================================================
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  optical_mat
    use od_constants,only : bohr2ang, H2eV
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
    character(len=100):: string, string2
    integer :: ik,is,ib,i,jb, ome_unit=6

    write(stdout,*) " Read a formatted ome file: "

    if(.not.allocated(optical_mat)) then
       write(stdout,*) "Allocating optical_mat"
       allocate(optical_mat(nbands,nbands,3,nkpoints,nspins))
    endif
    
    open(unit=ome_unit, form='formatted', file=trim(seedname)//".ome_fmt")

    ! Total number of elements of ome
    write(string,'(I0,"(x,",a,")")') 3*nbands*nbands, trim(format_precision)
    write(stdout,*) string

   ! write(string,'(a)') trim(format_precision)
    
    read(ome_unit,'('//trim(format_precision)//')') file_version
   
    read(ome_unit,'(a80)') omefile_header

   ! write(0,*) nkpoints, nspins, nbands
    
    do ik=1,nkpoints
       do is=1,nspins
          read(ome_unit,'('//trim(string)//')') (((optical_mat(ib,jb,i,ik,is),ib=1,nbands), &
               &jb=1,nbands),i=1,3)
       end do
    end do

    optical_mat=optical_mat*(bohr2ang*H2eV)

    close(unit=ome_unit)

     write(stdout,*) trim(seedname)//".ome_fmt"//"--> Formatted ome sucessfully read. "
    
  end subroutine read_ome_fmt
  
  !=========================================================================
  subroutine write_ome_fmt()
    !=========================================================================
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  optical_mat
    use od_constants,only : bohr2ang, H2eV
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
    character(len=100):: string
    integer :: ik,is,ib,i,jb, ome_unit=6


    optical_mat=optical_mat/(bohr2ang*H2eV)
    
    open(unit=ome_unit, form='formatted', file=trim(outseedname)//".ome_fmt")

    write(string,'(I0,"(x,",a,")")') 3*nbands*nbands, trim(format_precision)
    write(stdout,*) string

    write(stdout,'(a80)') omefile_header
    write(stdout,'(a80)') adjustl(omefile_header)
    
    write(ome_unit,'('//trim(format_precision)//')') file_version
    write(ome_unit,'(a80)') adjustl(omefile_header)
    
    do ik=1,nkpoints
       do is=1,nspins
          write(ome_unit,'('//trim(string)//')') (((optical_mat(ib,jb,i,ik,is),ib=1,nbands), &
               &jb=1,nbands),i=1,3)
       end do
    end do

    close(unit=ome_unit)

    write(stdout,*) " Sucesfully written a formatted ome file --> "//trim(outseedname)//".ome_fmt"
  end subroutine write_ome_fmt
  
  !=========================================================================
  subroutine read_ome_bin()
    !=========================================================================
    implicit none
    call elec_read_optical_mat()
    write(stdout,*) " "//trim(seedname)//".ome_bin"//"--> Unformatted ome sucessfully read. "
  end subroutine read_ome_bin
  
  !=========================================================================
  subroutine write_ome_bin()
    !=========================================================================
    use od_constants, only: dp, bohr2ang, H2eV
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  optical_mat
    use od_constants,only : bohr2ang, H2eV
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
    character(len=100):: string
    integer :: ik,is,ib,i,jb, ome_unit=6

    
    write(stdout,*) " Write a binary ome file"

    optical_mat=optical_mat/(bohr2ang*H2eV)
    
    open(unit=ome_unit, form='unformatted', file=trim(outseedname)//".ome_bin")

    write(stdout,*) "-> Omefile_version ", file_version
    write(ome_unit) file_version
    write(stdout,*) "-> Omefile_header ", trim(omefile_header)
    write(ome_unit) adjustl(omefile_header)
    
   ! write(0,*) nkpoints, nspins, nbands
    do ik=1,nkpoints
       do is=1,nspins
          write(ome_unit) (((optical_mat(ib,jb,i,ik,is),ib=1,nbands), &
               &jb=1,nbands),i=1,3)
       end do
    end do

     write(stdout,*) " Sucesfully written an unformatted ome file --> "//trim(outseedname)//".ome_bin"
  end subroutine write_ome_bin
  
  !=========================================================================
  ! D I A G O N A L  O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================
  
  !=========================================================================
  subroutine read_dome_fmt()
    !=========================================================================
       use od_constants, only: dp
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands, band_gradient
    use od_constants,only : bohr2ang, H2eV
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
    character(len=100):: string
    integer :: ik,is,ib,i,jb, dome_unit=6


    write(stdout,*) " Read a formatted dome file: "

    if(.not.allocated(band_gradient)) then
       write(0,*) " Allocating band_gradient"
       allocate(band_gradient(nbands,3,nkpoints,nspins))
    endif
    
    open(unit=dome_unit, form='formatted', file=trim(seedname)//".dome_fmt")
    
    write(string,'(i0,"(x,",a,")")') 3*nbands, trim(format_precision)
 
    read(dome_unit,'('//trim(format_precision)//')') file_version
   
    read(dome_unit,'(a80)') domefile_header
    
    do ik=1,nkpoints
       do is=1,nspins
          read(dome_unit,'('//trim(string)//')') ((band_gradient(ib,i,ik,is),ib=1,nbands), &
               &i=1,3)
       end do
    end do

    close(unit=dome_unit)

    
     write(stdout,*) " "//trim(seedname)//".dome_fmt"//"--> Formatted ome sucessfully read. "
  end subroutine read_dome_fmt
  
  !=========================================================================
  subroutine  write_dome_fmt()
    !=========================================================================
    use od_constants, only: dp
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  band_gradient
    use od_constants,only : bohr2ang, H2eV
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
    character(len=100):: string
    integer :: ik,is,ib,i,jb, dome_unit=6
    

    write(stdout,*) " Write a formatted ome file..."
    
    open(unit=dome_unit, form='formatted', file=trim(outseedname)//".dome_fmt")
    
    write(string,'(I0,"(x,",a,")")') 3*nbands, trim(format_precision)
    
    write(dome_unit,'('//trim(format_precision)//')') file_version
    write(dome_unit,'(a80)') adjustl(domefile_header)
    
    do ik=1,nkpoints
       do is=1,nspins
          write(dome_unit,'('//trim(string)//')') ((band_gradient(ib,i,ik,is),ib=1,nbands), &
            i=1,3)
       end do
    end do

    close(unit=dome_unit)

    write(stdout,*) " Sucesfully written a formatted dome file --> "//trim(outseedname)//".dome_fmt"
  end subroutine write_dome_fmt
  
  !=========================================================================
  subroutine read_dome_bin()
    !=========================================================================
    implicit none
    write(stdout,*) "Read a binary dome file"
    call elec_read_band_gradient()
  end subroutine read_dome_bin
  
  !=========================================================================
  subroutine  write_dome_bin()
    !=========================================================================
    
    use od_constants, only: dp
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  band_gradient
    implicit none
    
    real(dp):: file_version=1.0_dp          ! File version
 
    
    
    integer :: ik,is,ib,i,jb, dome_unit=6
    
    
    !  write(stdout,*) " Write a binary dome file"
    
    open(unit=dome_unit, form='unformatted', file=trim(outseedname)//".dome_bin")
    
    write(dome_unit) file_version
    write(dome_unit) adjustl(domefile_header)
    
    do ik=1,nkpoints
       do is=1,nspins
          write(dome_unit) ((band_gradient(ib,i,ik,is),ib=1,nbands), &
               &i=1,3)
       end do
    end do
  end subroutine write_dome_bin
  
  !=========================================================================
  ! P R O J E C T E D   D O S  
  !=========================================================================
  
  !=========================================================================
  subroutine  write_pdos_fmt()
    !=========================================================================
    use od_electronic, only : pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights 
    use od_cell, only : nkpoints, kpoint_r
    use od_io,only : stdout, seedname, io_file_unit
    implicit none

    integer :: ik,is,ib
    integer :: pdos_in_unit
    character(len=80) :: string
    real(dp) :: file_version=1.0_dp
    
    !-------------------------------------------------------------------------!
    ! W R I T E   T H E   D A T A   H E A D E R

   
    pdos_in_unit=io_file_unit()

    open(unit=pdos_in_unit,file=trim(seedname)//".pdos_fmt",form='formatted')
    write(pdos_in_unit,'('//trim(format_precision)//')') file_version
    write(pdos_in_unit,'(a80)') trim(pdosfile_header)
    
    write(pdos_in_unit,'(i6)') pdos_mwab%nkpoints
    write(pdos_in_unit,'(i6)') pdos_mwab%nspins
    write(pdos_in_unit,'(i6)') pdos_mwab%norbitals
    write(pdos_in_unit,'(i6)') pdos_mwab%nbands
   
 
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nkpoints= ",pdos_mwab%nkpoints
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nspins= ",pdos_mwab%nspins
    write(stdout,'(i6)') "DEBUG: pdos_mwab%norbitals= ",pdos_mwab%norbitals
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nbands= ",pdos_mwab%nbands


    ! These should all be allocated!
    !allocate(pdos_orbital%species_no(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    !allocate(pdos_orbital%rank_in_species(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    !allocate(pdos_orbital%am_channel(pdos_mwab%norbitals),stat=ierr)
    !if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")

    write(string,'(i7,"(x,"a")")') pdos_mwab%norbitals, trim(format_precision)
   
    
    write(pdos_in_unit,'('//trim(string)//')') pdos_orbital%species_no(1:pdos_mwab%norbitals)
    write(pdos_in_unit,'('//trim(string)//')') pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    write(pdos_in_unit,'('//trim(string)//')') pdos_orbital%am_channel(1:pdos_mwab%norbitals)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A

    ! These should already be alloacted
    ! allocate(nbands_occ(1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate nbands_occ"
    ! allocate(pdos_weights(1:pdos_mwab%norbitals,1:pdos_mwab%nbands, &
    !     1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate pdos_weights"
 
    do ik=1,nkpoints
       ! The kpoint number, followed by the kpoint-vector
       write(pdos_in_unit,'(i6,3'//trim(format_precision)//')') ik, kpoint_r(:,ik)
       do is=1, pdos_mwab%nspins
          write(pdos_in_unit,'(i6)') is ! this is the spin number
          write(pdos_in_unit,'(i6)') nbands_occ(ik,is)
          do ib=1,nbands_occ(ik,is)
             
       !      write(stdout,*) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
       !      write(stdout,*) " DEBUG:", ib, ik, is
       !      write(stdout,*) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "
             
             write(pdos_in_unit,'('//trim(string)//')') pdos_weights(1:pdos_mwab%norbitals,ib,ik,is)
          enddo
       enddo
    enddo
           
    close(pdos_in_unit)
   
  end subroutine write_pdos_fmt
  
  !=========================================================================
  subroutine  read_pdos_fmt()
    !=========================================================================
    use od_electronic, only : pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights 
    use od_cell, only : nkpoints, kpoint_r
    use od_io,only : stdout, seedname, io_file_unit
    implicit none

    integer :: ik,is,ib, idummy, ierr
    integer :: pdos_in_unit
    character(len=80) :: string
    real(dp) :: file_version
    
    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R


    file_version=1.0_dp
    pdos_in_unit=io_file_unit()

    open(unit=pdos_in_unit,file=trim(seedname)//".pdos_fmt",form='formatted')
    read(pdos_in_unit,'('//trim(format_precision)//')') file_version
    read(pdos_in_unit,'(a80)') pdosfile_header
    
    read(pdos_in_unit,'(i6)') pdos_mwab%nkpoints
    read(pdos_in_unit,'(i6)') pdos_mwab%nspins
    read(pdos_in_unit,'(i6)') pdos_mwab%norbitals
    read(pdos_in_unit,'(i6)') pdos_mwab%nbands
   
 
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nkpoints= ",pdos_mwab%nkpoints
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nspins= ",pdos_mwab%nspins
    write(stdout,'(i6)') "DEBUG: pdos_mwab%norbitals= ",pdos_mwab%norbitals
    write(stdout,'(i6)') "DEBUG: pdos_mwab%nbands= ",pdos_mwab%nbands


    write(string,'(i7,"(x,",a,")")') pdos_mwab%norbitals, trim(format_precision)
    
    ! These should all be allocated!
    allocate(pdos_orbital%species_no(pdos_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    allocate(pdos_orbital%rank_in_species(pdos_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")
    allocate(pdos_orbital%am_channel(pdos_mwab%norbitals),stat=ierr)
    if(ierr/=0) call io_error(" Error : cannot allocate pdos_orbital")


   
    
    read(pdos_in_unit,'('//trim(string)//')') pdos_orbital%species_no(1:pdos_mwab%norbitals)
    read(pdos_in_unit,'('//trim(string)//')') pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    read(pdos_in_unit,'('//trim(string)//')') pdos_orbital%am_channel(1:pdos_mwab%norbitals)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A

    ! These should already be alloacted
    ! allocate(nbands_occ(1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate nbands_occ"
    ! allocate(pdos_weights(1:pdos_mwab%norbitals,1:pdos_mwab%nbands, &
    !     1:num_kpoints_on_node(my_node_id),1:pdos_mwab%nspins),stat=ierr)
    ! if(ierr/=0) stop " Error : cannot allocate pdos_weights"
 
    do ik=1,nkpoints
       ! The kpoint number, followed by the kpoint-vector
       read(pdos_in_unit,'(i6,3'//trim(format_precision)//')') idummy, kpoint_r(:,ik)
       do is=1, pdos_mwab%nspins
          read(pdos_in_unit,'(i6)') idummy ! this is the spin number
          read(pdos_in_unit,'(i6)') nbands_occ(ik,is)
          do ib=1,nbands_occ(ik,is)
             
       !      write(stdout,*) " ***** F U L L _ D E B U G _ P D O S _ W E I G H T S ***** "
       !      write(stdout,*) " DEBUG:", ib, ik, is
       !      write(stdout,*) "   **** ***** *****  ***** ***** *****  ***** ***** *****  "
             
             read(pdos_in_unit,'('//trim(string)//')') pdos_weights(1:pdos_mwab%norbitals,ib,ik,is)
          enddo
       enddo
    enddo
           
    close(pdos_in_unit)

  end subroutine read_pdos_fmt
  
  !=========================================================================
  subroutine  read_pdos_bin()
    !=========================================================================
    implicit none
    write(stdout,*) "Read a binary pdos file"
    call elec_pdos_read()
  end subroutine read_pdos_bin
  
  !=========================================================================
  subroutine  write_pdos_bin()
    !=========================================================================
    use od_constants, only : dp
    use od_electronic, only : pdos_mwab, nbands_occ, pdos_orbital, nspins, pdos_weights 
    use od_cell, only : nkpoints, kpoint_r
    use od_io, only : seedname
    
    implicit none
    
  !  integer,parameter:: num_popn_orb=1          ! Number of LCAO projectors
   ! integer,parameter:: max_eigenv=1       ! Number of bands included in matrix elements
    real(dp):: file_version=1.0_dp  ! File version
   ! integer:: species(1:num_popn_orb)! Atomic species associated with each projector
   ! integer:: ion(1:num_popn_orb)    ! Ion associated with each projector
   ! integer:: am_channel(1:num_popn_orb)     ! Angular momentum channel
   ! integer:: num_eigenvalues(1:nspins)   ! Number of eigenvalues per spin channel
   ! real(dp):: kpoint_positions(1:nkpoints,1:3) ! k_x, k_y, k_z in fractions of BZ
   ! real(dp):: pdos_weights(1:num_popn_orb,max_eigenv,num_kpoints,num_spins)!Matrix elements
    !character(len=80):: file_header ! File header comment
    
    integer :: ik,is,ib, pdos_file=6
    
    open(unit=pdos_file, form='unformatted', file="pdos.out")
    
    write(stdout,*) " Write a binary pdos file"
    
    write(pdos_file) file_version
    write(pdos_file) trim(pdosfile_header)
    write(pdos_file) pdos_mwab%nkpoints
    write(pdos_file) pdos_mwab%nspins
    write(pdos_file) pdos_mwab%norbitals
    write(pdos_file) pdos_mwab%nbands
    write(pdos_file) pdos_orbital%species_no(1:pdos_mwab%norbitals)
    write(pdos_file) pdos_orbital%rank_in_species(1:pdos_mwab%norbitals)
    write(pdos_file) pdos_orbital%am_channel(1:pdos_mwab%norbitals)
    
    do ik=1,pdos_mwab%nkpoints
       write(pdos_file) ik, kpoint_r(ik,:)
       do is = 1,pdos_mwab%nspins
          write(pdos_file) is
          write(pdos_file) nbands_occ(ik,is)
          do ib = 1,nbands_occ(ik,is)
             write(pdos_file) real(pdos_weights(pdos_mwab%norbitals,ib,ik,is))
          end do
       end do
    end do
    
    
  end subroutine write_pdos_bin
  
  !=========================================================================
  ! E L N E S   M A T R I X   E L E M E N T S
  !=========================================================================
  
  !=========================================================================
  subroutine  read_elnes_fmt()
    !=========================================================================
    implicit none
    write(stdout,*) "Read a formatted elnes file"
    write(stdout,*) "Not implemented"
  end subroutine read_elnes_fmt
  
  !=========================================================================
  subroutine  write_elnes_fmt()
    !=========================================================================
    implicit none
    write(stdout,*) "Write a formatted elnes file"
    write(stdout,*) "Not implemented"
  end subroutine write_elnes_fmt

    !=========================================================================
  subroutine  read_elnes_bin()
    !=========================================================================
    implicit none
    write(stdout,*) "Read a binary elnes file"
    call  elec_read_elnes_mat()
  end subroutine read_elnes_bin
  
  !=========================================================================
  subroutine  write_elnes_bin()
    !=========================================================================
    implicit none
    
    integer,parameter:: dp=selected_real_kind(15,300)
    integer,parameter:: tot_core_projectors=1  ! Total number of core states included in matrix elements
    integer,parameter:: max_eigenvalues=1      ! Number of bands included in matrix elements
    integer,parameter:: num_kpoints=1          ! Number of k-points
    integer,parameter:: num_spins=1            ! Number of spins
    integer,parameter:: num_eigenvalues(1:num_spins)=1    ! Number of eigenvalues
    integer:: species(1:tot_core_projectors) ! Atomic species associated with each projector
    integer:: ion(1:tot_core_projectors)     ! Ion associated with each projector
    integer:: n(1:tot_core_projectors)       ! Principal quantum number associated with each projector
    integer:: lm(1:tot_core_projectors)      ! Angular momentum quantum number associated with each projector
    real(dp):: elnes_mat(tot_core_projectors,max_eigenvalues,3,num_kpoints,num_spins) ! matrix elements
    
    
    integer :: nk,ns,nb,i,jb, elnes_file=6, orb,indx
    open(unit=elnes_file, form='unformatted', file="elnes.out")
    
    write(stdout,*) " Write a binary elnes file"

    !!write !! Some headers here?
    
    write(elnes_file) tot_core_projectors
    write(elnes_file) max_eigenvalues
    write(elnes_file) num_kpoints
    write(elnes_file) num_spins
    write(elnes_file) species(1:tot_core_projectors)
    write(elnes_file) ion(1:tot_core_projectors)
    write(elnes_file) n(1:tot_core_projectors)
    write(elnes_file) lm(1:tot_core_projectors)
    
    do nk = 1,num_kpoints
       do ns = 1, num_spins
          write(elnes_file) (((elnes_mat(orb,nb,indx,nk,ns),orb=1,&
               &tot_core_projectors),nb=1,num_eigenvalues(ns)),indx=1,3)
       end do
    end do
    
    
  end subroutine write_elnes_bin
  

  !=========================================================================
  subroutine slice_an_ome()
    !=========================================================================
    use od_constants, only: dp
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  band_gradient, optical_mat
    implicit none
    
    integer :: loop
    
    write(stdout,*) " Slicing an ome into a dome."
    
    if(.not.allocated(band_gradient)) then
       write(0,*) " Allocating band_gradient"
       allocate(band_gradient(nbands,3,nkpoints,nspins))
    endif
    
    do loop=1,nbands
       band_gradient(loop,:,:,:)=real(optical_mat(loop,loop,:,:,:),dp)
    end do
    
  end subroutine slice_an_ome

  !=========================================================================
  subroutine pad_an_ome()
    !=========================================================================
    use od_constants, only: dp
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_electronic, only : nspins, nbands,  band_gradient, optical_mat
    implicit none
    
    integer :: loop
    
    write(stdout,*) " Padding a dome into a ome."
    
    if(.not.allocated(optical_mat)) then
       write(0,*) " Allocating optical_mat"
       allocate(optical_mat(nbands,nbands,3,nkpoints,nspins))
    endif

    ! Pad with zeros.
    optical_mat=0.0_dp
    
    do loop=1,nbands
       optical_mat(loop,loop,:,:,:)=band_gradient(loop,:,:,:)
    end do
    
  end subroutine pad_an_ome

  !=========================================================================
  subroutine  report_arraysize()
    !=========================================================================
    use od_electronic, only:  nspins, nbands
    use od_cell, only: nkpoints
    use od_io, only: stdout
    implicit none

    write(stdout,*)
    write(stdout,'(a40,i5)') " Number of kpoints : ", nkpoints
    write(stdout,'(a40,i5)') " Number of bands : ", nbands
    write(stdout,'(a40,i5)') " Number of spins : ", nspins
    write(stdout,*)
    
  end subroutine report_arraysize
  
  end module od_conv
          


program od2od
  !! Program to convert checkpoint files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use od_constants, only: dp
  use od_io, only: io_file_unit, stdout, stderr, io_error, seedname, io_time,io_date
  use od_conv
  use od_comms, only: num_nodes, comms_setup, comms_end
  implicit none

  logical :: file_found
  logical ::   ome_conv, dome_conv, pdos_conv, elnes_conv ! Flags to stop people trying
  ! to read in a pdos and write out an elnes.  That's not going to end well.
  integer :: file_unit
  real(kind=dp)    :: time0,time1       ! Varaibles for timing
  character(len=9) :: stat,pos          ! Status and position of .odo file
  character(len=9) :: ctime             ! Temp. time string
  character(len=11):: cdate             ! Temp. date string

  time0=io_time()

  iprint=4
  
  call comms_setup

 
  call conv_get_seedname
  stderr=io_file_unit()
  open(unit=stderr,file=trim(seedname)//'.opt_err')
  call io_date(cdate,ctime)
  write(stderr,*)  'od2od: Execution started on ',cdate,' at ',ctime
  
  stdout = io_file_unit()
  open (unit=stdout, file=trim(seedname)//'.log')
  !-------------------------------------------------------------------------!

  write(stdout,*)
  write(stdout,*) " +=======================================================================+ "
  write(stdout,*) " |                                                                       | "
  write(stdout,*) " |                      OO   DDD   222    OO   DDD                       | "
  write(stdout,*) " |                     O  O  D  D     2  O  O  D  D                      | "
  write(stdout,*) " |                     O  O  D  D   22   O  O  D  D                      | "
  write(stdout,*) " |                     O  O  D  D  2     O  O  D  D                      | "
  write(stdout,*) " |                      OO   DDD   2222   OO   DDD                       | "
  write(stdout,*) " |                                                                       | "
  write(stdout,*) " |                For doing the odd thing to OptaDOS files               | "
  write(stdout,*) " |                                                                       | "
  write(stdout,*) " |                   OptaDOS Developers Group 2019 (C)                   | "                                  
  write(stdout,*) " |                          (But blame Andrew)                           | "
  write(stdout,*) " |                                                                       | "
  write(stdout,*) " +=======================================================================+ "

  if (num_nodes /= 1) then
    call io_error('od2od can only be used in serial...')
  endif

  write(stdout,*) 
  write(stdout,'(a40,i5)') "Number of nodes : ",num_nodes
  write(stdout,'(a40,a)') "Convert from : ", trim(infile)
  write(stdout,'(a40,a)') "Convert to : ", trim(outfile)
  write(stdout,'(a40,a)') "Seedname : ", trim(seedname)
  write(stdout,'(a40,a)') "Output Seedname : ", trim(outseedname)



  ome_conv=.false.
  dome_conv=.false.
  pdos_conv=.false.
  elnes_conv=.false.
  

  ! Main case to decide what file format to read in.
  read_input: select case(trim(infile))
  case("ome_fmt")
     ome_conv=.true.
     call elec_read_band_energy()
     call report_arraysize()
     call read_ome_fmt()
  case("ome_bin")
     ome_conv=.true.
     call elec_read_band_energy()
     call report_arraysize()
     call read_ome_bin()
  case ("dome_fmt")
     dome_conv=.true.
     call elec_read_band_energy()
     call report_arraysize()
     call read_dome_fmt()
  case("dome_bin")
     dome_conv=.true.
     call elec_read_band_energy()
     call report_arraysize()
     call read_dome_bin()
     write(stdout,*) domefile_header
  case ("pdos_fmt")
     pdos_conv=.true.
     call elec_read_band_energy()
     call report_arraysize()
     call read_pdos_fmt()
  case("pdos_bin")
      pdos_conv=.true.
      call elec_read_band_energy()
      call report_arraysize()
     call read_pdos_bin()
  case ("elnes_fmt")
      elnes_conv=.true.
      call elec_read_band_energy()
      call report_arraysize()
     call read_elnes_fmt()
  case("elnes_bin")
      elnes_conv=.true.
      call elec_read_band_energy()
      call report_arraysize()
     call read_elnes_bin()
  case default
     call io_error('Unknown Input File format speccified')
  end select read_input

  ! Main case to decide what file format to write.
  write_output: select case(trim(outfile))
  case("ome_fmt")
     if(.not. (dome_conv.or.ome_conv)) call io_error(' Input format '//trim(infile)//'not compatible with output format '&
          &//trim(outfile))
     if(dome_conv) call pad_an_ome() 
     call write_ome_fmt() 
  case("ome_bin")
     if(.not.  (dome_conv.or.ome_conv)) call io_error(' Input format '//trim(infile)//'not compatible with output format '&
          &//trim(outfile))
     if(dome_conv) call pad_an_ome() 
     call write_ome_bin() 
  case ("dome_fmt")
     if(.not. (dome_conv.or.ome_conv)) call io_error(' Input format '//trim(infile)//&
          &'not compatible with output format '//trim(outfile))
     if(ome_conv) call slice_an_ome()    
     call write_dome_fmt()
  case("dome_bin")
     if(.not. (dome_conv.or.ome_conv)) call io_error(' Input format '//trim(infile)//&
          &'not compatible with output format '//trim(outfile))
     if(ome_conv) call slice_an_ome()
     call write_dome_bin() 
  case ("pdos_fmt")
      if(.not. pdos_conv) call io_error(' Input format '//trim(infile)//'not compatible with output format '//trim(outfile))
     call write_pdos_fmt() 
  case("pdos_bin")
     if(.not. pdos_conv) call io_error(' Input format '//trim(infile)//'not compatible with output format '//trim(outfile))
     call write_pdos_bin()  
  case ("elnes_fmt")
     if(.not. elnes_conv) call io_error(' Input format '//trim(infile)//'not compatible with output format '//trim(outfile))
     call write_elnes_fmt()  
  case("elnes_bin")
     if(.not. elnes_conv) call io_error(' Input format '//trim(infile)//'not compatible with output format '//trim(outfile))
     call write_elnes_bin()  
  case default
     call io_error('Unknown Output File format speccified')
  end select write_output


  call io_date(cdate,ctime)
  write(stdout,*)
  time1=io_time()
  
  write(stdout,'(1x,a40,f11.3,a)') 'Total runtime ',time1-time0,' (sec)'
  write(stdout,*)
  write(stdout,*) 'od2od: Execution complete on ',cdate,' at ',ctime
  write(stdout,*)
  
  close(stdout)
  close(stderr, status='delete')
  
  
  call comms_end

 
!  close(unit=stdout,status='delete')


end program od2od
