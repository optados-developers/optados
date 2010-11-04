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
program optados
  use constants, only : dp
  use od_io, only : io_get_seedname, stdout, io_time, io_date, io_file_unit, seedname
  use od_parameters
  use dos ! only : ??
  implicit none

  logical :: odo_found
  character(len=9) :: stat,pos,cdate,ctime
  real(kind=dp) :: time0,time1,time2


  time0=io_time()
  call io_get_seedname()

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.opt_err')
  call io_date(cdate,ctime)
  write(stdout,*)  'OptaDos: Execution started on ',cdate,' at ',ctime
  call param_read()
  close(stdout,status='delete')

  inquire(file=trim(seedname)//'.odo',exist=odo_found)
  if (odo_found) then
     stat='old'
  else
     stat='replace'
  endif
  pos='append'

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.odo',status=trim(stat),position=trim(pos))
!  call param_write_header()
  call param_write()

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'



 call param_dealloc

endprogram
