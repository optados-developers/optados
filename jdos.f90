!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
module od_jdos
  use od_constants, only : dp
  implicit none
  
    !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S 
  public :: jdos_calculate
  
  contains

!===============================================================================  
  subroutine jdos_calculate
  !===============================================================================
   use od_jdos_utils, only : jdos_utils_calculate, write_jdos,E,jdos_fixed &
   &,jdos_adaptive, jdos_linear
   use od_parameters, only : fixed, adaptive, linear
   use od_io, only : io_time,stdout
   use od_comms, only : on_root
   implicit none
   
   real(dp) :: time0,time1
   call jdos_utils_calculate

!-------------------------------------------------------------------------------
! W R I T E   O U T   D O S  

time0=io_time()
! Otherwise we have written to wdos and dos, so they can be called 
! by whatever.
  if(on_root) then
    if(fixed)    call write_jdos(E, jdos_fixed,  "fixed")
    if(adaptive) call write_jdos(E, jdos_adaptive, "adaptive")
    if(linear)   call write_jdos(E, jdos_linear,  "linear")
    !if(quad)    call write_jdos(E, dos_quad, intdos_quad, "quad")
   endif
time1=io_time()
write(stdout,'(1x,a40,f11.3,a)') 'Time to write jdos to disk ',time1-time0,' (sec)'

!-------------------------------------------------------------------------------

  end subroutine jdos_calculate
  !===============================================================================
endmodule od_jdos