!-*- mode: F90 -*-!
module  od_conv
  use od_constants, only: dp
  use od_electronic, only: elec_read_optical_mat, elec_read_band_gradient,  elec_read_elnes_mat, elec_pdos_read 
  use od_io, only:  stdout, io_error, seedname
  implicit none

  logical, save :: format2unformatted

 contains
 !=========================================================================
   subroutine print_usage()
 !=========================================================================
    !! Writes the usage of the program to stdout
    write (stdout, '(A)') " When I know what's happening, you'll be the first to know...."
  end subroutine print_usage 

  !=========================================================================
  subroutine conv_get_seedname
    !=========================================================================
    !! Set the seedname from the command line
    implicit none
    
    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 1) then
      seedname = 'optados'
    elseif (num_arg == 2) then
      call get_command_argument(2, seedname)
    else
      call print_usage
      call io_error('Wrong command line arguments, see logfile for usage')
    end if

    call get_command_argument(1, ctemp)
    if (index(ctemp, '-f2u') > 0) then
      format2unformatted = .false.
    elseif (index(ctemp, '-u2f') > 0) then
       format2unformatted = .true.
    else
      write (stdout, '(A)') 'Wrong command line action: '//trim(ctemp)
      call print_usage
      call io_error('Wrong command line arguments, see logfile for usage')
    end if

  end subroutine conv_get_seedname 

  !=========================================================================
  ! O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================

  
  !=========================================================================
  subroutine read_ome_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Read in a formatted ome file"
    write(*,*) "Not implemented"
  end subroutine read_ome_fmt
  
  !=========================================================================
  subroutine write_ome_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Write a formatted ome file"
    write(*,*) "Not implemented"
  end subroutine write_ome_fmt
  
  !=========================================================================
  subroutine read_ome_bin()
    !=========================================================================
    implicit none
    write(*,*) " Read a binary ome file..."
    call elec_read_optical_mat()
  end subroutine read_ome_bin
  
  !=========================================================================
  subroutine write_ome_bin()
    !=========================================================================
    implicit none
    write(*,*) "Write a binary ome file"
    write(*,*) "We wrote a snippet of code to do this"
  end subroutine write_ome_bin
  
  !=========================================================================
  ! D I A G O N A L  O P T I C A L   M A T R I X   E L E M E N T S
  !=========================================================================
  
  !=========================================================================
  subroutine read_dome_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Read a binary dome file"
    write(*,*) "Not implemented"
  end subroutine read_dome_fmt
  
  !=========================================================================
  subroutine  write_dome_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Write a binary dome file"
    write(*,*) "Not implemented"
  end subroutine write_dome_fmt
  
  !=========================================================================
  subroutine read_dome_bin()
    !=========================================================================
    implicit none
    write(*,*) "Read a binary dome file"
    call elec_read_band_gradient()
  end subroutine read_dome_bin
  
  !=========================================================================
  subroutine  write_dome_bin()
    !=========================================================================
    implicit none
    write(*,*) "Write a binary dome file"
    write(*,*) "We wrote a snippet of code to do this"
    write(*,*) "Not implemented"
  end subroutine write_dome_bin
  
  !=========================================================================
  ! P R O J E C T E D   D O S  
  !=========================================================================
  
  !=========================================================================
  subroutine  read_pdos_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Read a formatted dome file"
    write(*,*) "Not implemented"
  end subroutine read_pdos_fmt
  
  !=========================================================================
  subroutine  write_pdos_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Write a formatted dome file"
    write(*,*) "Not implemented"
  end subroutine write_pdos_fmt
  
  !=========================================================================
  subroutine  read_pdos_bin()
    !=========================================================================
    implicit none
    write(*,*) "Read a binary pdos file"
    call elec_pdos_read()
  end subroutine read_pdos_bin
                    
  !=========================================================================
  subroutine  write_pdos_bin()
    !=========================================================================
        implicit none
        write(*,*) "Write a binary pdos file"
        write(*,*) "We wrote a snippet of code to do this"
    write(*,*) "Not implemented"
  end subroutine write_pdos_bin

    !=========================================================================
  ! E L N E S   M A T R I X   E L E M E N T S
  !=========================================================================
  
  !=========================================================================
  subroutine  read_elnes_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Read a formatted elnes file"
    write(*,*) "Not implemented"
  end subroutine read_elnes_fmt
  
  !=========================================================================
  subroutine  write_elnes_fmt()
    !=========================================================================
    implicit none
    write(*,*) "Write a formatted elnes file"
    write(*,*) "Not implemented"
  end subroutine write_elnes_fmt

    !=========================================================================
  subroutine  read_elnes_bin()
    !=========================================================================
    implicit none
    write(*,*) "Read a binary elnes file"
    call  elec_read_elnes_mat()
  end subroutine read_elnes_bin
  
  !=========================================================================
  subroutine  write_elnes_bin()
    !=========================================================================
    implicit none
    write(*,*) "Write a binary elnes file"
    write(*,*) "We wrote a snippet of code to do this"
    write(*,*) "Not implemented"
  end subroutine write_elnes_bin
  

  

  
  end module od_conv
          


program od_fmt2bin
  !! Program to convert checkpoint files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use od_constants, only: dp
  use od_io, only: io_file_unit, stdout, io_error, seedname
  use od_conv
  use od_comms, only: num_nodes, comms_setup, comms_end
  implicit none

  ! Export mode:
  !  TRUE:  create formatted .chk.fmt from unformatted .chk ('-export')
  !  FALSE: create unformatted .chk from formatted .chk.fmt ('-import')
  logical :: file_found
  integer :: file_unit

  call comms_setup

  stdout = io_file_unit()
  open (unit=stdout, file='od_fmt2bin.log')

  if (num_nodes /= 1) then
    call io_error('of_fmt2bin can only be used in serial...')
  endif

  call conv_get_seedname

 ! if (export_flag .eqv. .true.) then
 !   call conv_read_chkpt()
 !   call conv_write_chkpt_fmt()
 ! else
 !   call conv_read_chkpt_fmt()
 !   call conv_write_chkpt()
  ! end if

  if ( format2unformatted .eqv. .true.) then
     write(*,*)  " Going to turn formatted files into binaries."
     call read_ome_bin()
     call write_ome_fmt()
     call read_dome_bin()
     call write_dome_fmt()
     call read_pdos_bin()
     call write_pdos_fmt()
     call read_elnes_bin()
     call write_elnes_fmt()
  elseif ( format2unformatted .eqv. .false. ) then
     write(*,*) " Going to turn binary files into formatted ones."
     call read_ome_fmt()
     call write_ome_bin()
     call read_dome_fmt()
     call write_dome_bin()
     call read_pdos_fmt()
     call write_pdos_bin()
     call read_elnes_fmt()
     call write_elnes_bin()
  endif

!  close(unit=stdout,status='delete')
  close (unit=stdout)

  call comms_end

end program od_fmt2bin
