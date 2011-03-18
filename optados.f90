program optados
!=========================================================================!
!                             O P T A D O S                               !
!                      OPTics And Density Of States                       !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below.                                                        !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris, Rebecca Nicholls, Chris Pickard               ! 
!             and Jonathan Yates      2010                                !
!=========================================================================!
  use od_comms, only  : comms_setup,on_root,comms_end,my_node_id
  use od_constants, only : dp
  use od_io, only     : io_get_seedname, io_time, io_date, io_file_unit,&! Functions
   & stdout, stderr, seedname                                            ! Variables
  use od_parameters
  use od_cell
  use od_electronic,only : elec_read_band_energy,elec_report_parameters
  use od_dos,  only : dos_calculate
  use od_jdos,  only: jdos_calculate
  use od_pdos, only : dos_partial, pdos_write
!  use od_pdos, only : pdos_write, pdos_weights, dos_partial
  implicit none

  real(kind=dp)    :: time0,time1,time2 ! Varaibles for timing
  logical          :: odo_found         ! Ouptut file exists?
  character(len=9) :: stat,pos          ! Status and position of .odo file
  character(len=9) :: ctime             ! Temp. time string
  character(len=11):: cdate             ! Temp. date string

  time0=io_time()

  call comms_setup

  if (on_root) then
     call io_get_seedname()
     !-------------------------------------------------------------------------!
     ! O R G A N I S E   T H E   E R R O R   F I L E 
     stderr=io_file_unit()
     open(unit=stderr,file=trim(seedname)//'.opt_err')
     call io_date(cdate,ctime)
     write(stderr,*)  'OptaDOS: Execution started on ',cdate,' at ',ctime
     !-------------------------------------------------------------------------!
     
     
     !-------------------------------------------------------------------------!
     ! O R G A N I S E   T H E   O U T P U T   F I L E  A N D 
     ! R E A D   A N D   W R I T E   U S E R   P A R A M E T E R S   
     call param_read()
     inquire(file=trim(seedname)//'.odo',exist=odo_found)
     if (odo_found) then
        stat='old'
     else
        stat='replace'
     endif
     pos='append'

     stdout=io_file_unit()
     open(unit=stdout,file=trim(seedname)//'.odo',status=trim(stat),position=trim(pos))
     write(stdout,*)  'OptaDOS: Execution started on ',cdate,' at ',ctime
     call param_write_header()
     call param_write()
     time1=io_time()
     write(stdout,*)
     write(stdout,'(1x,a40,f11.3,a)') 'Time to read parameters ',time1-time0,' (sec)'
     !-------------------------------------------------------------------------!

     call elec_read_band_energy
     call cell_calc_lattice
     call cell_report_parameters
     call elec_report_parameters
    
  end if
  ! now send the data from the parameter file to each node
 
  call param_dist
  call cell_dist
 
!-------------------------------------------------------------------------!
! C A L L   P D O S   R O U T I N E S
  if(pdos) then
    time0=io_time()
!   call cell_pdos_read 
!   call pdos_merge
    time1=io_time()
    write(stdout,'(1x,a40,f11.3,a)') 'Time to set up Partial DOS ',time1-time0,' (sec)'

    time0=io_time()
!    call dos_calculate(matrix_weights=pdos_weights, weighted_dos=dos_partial)
    call pdos_write

    time1=io_time()
    if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Partical DOS ',time1-time0,' (sec)'
  endif
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! C A L L   C O R E   R O U T I N E S
  if(core) then
    time0=io_time()
    !call core_calculate
    time1=io_time()
    if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Core Level Spec. (Total) ',time1-time0,' (sec)'
  endif
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! C A L L   D O S  R O U T I N E S
  if(dos) then
    time0=io_time()
    call dos_calculate
    time1=io_time()
    if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate DOS (Total) ',time1-time0,' (sec)'
  endif
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! C A L L   O P T I C S   R O U T I N E S
  if(optics) then
    time0=io_time()
    !call optics_calculate
    time1=io_time()
      if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Optical spec. DOS (Total) ',time1-time0,' (sec)'
  endif
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! C A L L   J D O S   R O U T I N E S
  if(jdos) then
    time0=io_time()
    call jdos_calculate
    time1=io_time()
      if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Joint DOS (Total) ',time1-time0,' (sec)'
  endif
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! F I N A L I S E 
   call param_dealloc

   if(on_root) then
      call io_date(cdate,ctime)
      write(stdout,*)
      write(stdout,*) 'OptaDOS: Execution complete on ',cdate,' at ',ctime
      
      close(stdout)
      close(stderr, status='delete')
   end if

   call comms_end

!-------------------------------------------------------------------------!
end program optados
