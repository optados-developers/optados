!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
module od_dos
  use od_constants, only : dp
  use od_dos_utils, only : dos_utils_calculate
  implicit none
  
    !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S 
  public :: dos_calculate
  
  contains
  
  subroutine dos_calculate
    !===============================================================================  
    ! Main routine in dos module, drives the calculation of density of states for
    ! both task : dos and also if it is required elsewhere.
    !------------------------------------------------------------------------------- 
    ! Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
    !            weighted_dos   (out)(opt) : Output DOS weigthed by matrix_weights
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
    ! intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
    ! efermi_linear, delta_bins, calc_weighted_dos
    !-------------------------------------------------------------------------------
    ! Modules Used: see below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: One of linear, adaptive or fixed must be .true.
    !-------------------------------------------------------------------------------
    ! Known Worries: (1) If more than one of linear, adaptive or fixed are set it 
    ! uses the most complicated method.
    ! (2) It should be possible to pass optioinal arguments to sub programs as
    ! optional argumnets without checking whether they are there or not. g95 will 
    ! allow this behaviour. gfotran will not.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !=============================================================================== 
    use od_io,        only : stdout,io_time,io_error
    use od_dos_utils, only : E, dos_fixed, intdos_fixed, dos_adaptive, &
    &intdos_adaptive, dos_linear, intdos_linear
    use od_parameters,only : fixed,adaptive,linear
    use od_comms,     only : on_root
       
       
    real(dp) :: time0, time1
       
    call dos_utils_calculate   ! Will return if this has already been done.
       
    ! W R I T E   O U T   D O S  
       time0=io_time()
       ! Otherwise we have written to wdos and dos, so they can be called 
       ! by whatever.
       if(on_root) then
          if(fixed)    call write_dos(E, dos_fixed, intdos_fixed, "fixed")
          if(adaptive) call write_dos(E, dos_adaptive, intdos_adaptive, "adaptive")
          if(linear)   call write_dos(E, dos_linear,  intdos_linear, "linear")
          !if(quad)    call write_dos(E, dos_quad, intdos_quad, "quad")
       endif
       time1=io_time()
       if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to write dos to disk ',time1-time0,' (sec)'
  
    !-------------------------------------------------------------------------------
    
 end subroutine dos_calculate
 
  !=============================================================================== 
  subroutine write_dos(E,dos,intdos,dos_name)
    !=============================================================================== 
    ! This routine receives an energy scale, a density of states and a file name
    ! and writes out the DOS to disk
    !------------------------------------------------------------------------------- 
    ! Arguments: E       (in) : The energy scale
    !            dos     (in) : The density of states
    !            intdos  (in) : The integrated DOS
    !            dos_name(in) : Name of the output file
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: None
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !=============================================================================== 
    use od_electronic, only : nspins
    use od_parameters, only : nbins, dos_per_volume
    use od_io,         only : seedname, io_file_unit,io_date,io_error

    implicit none
    real(dp), intent(in) :: E(nbins)
    real(dp), intent(in) :: dos(nbins,nspins)
    real(dp), intent(in) :: intdos(nbins,nspins)
    character(len=*), intent(in) :: dos_name
    integer :: i, dos_file, ierr
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: dos_units, intdos_units


    dos_file=io_file_unit()
    open(unit=dos_file,file=trim(seedname)//'.'//trim(dos_name)//'.dat',iostat=ierr)
    if(ierr.ne.0) call io_error(" ERROR: Cannot open output file in dos: write_dos")

    dos_units="(electrons per eV)" ; intdos_units="(electrons)"
    if(dos_per_volume) then
       dos_units="(electrons per eV/A^3)" 
       intdos_units="(electrons per A^3)"
    endif

    write(dos_file, *) "##############################################################################"
    write(dos_file,*) "#"
    write(dos_file, *) "#                  O p t a D O S   o u t p u t   f i l e "  
    write(dos_file, '(1x,a1)') "#"
    write(dos_file,*) "#    Denisty of States using ", trim(dos_name), " broadening"
    call io_date(cdate,ctime)
    write(dos_file,*)  '#  Generated on ',cdate,' at ',ctime
    write(dos_file,*) "# Column        Data"
    write(dos_file,*) "#    1        Energy (eV)"
    if(nspins>1) then
       write(dos_file,*) "#    1        Up-spin DOS ", trim(dos_units)
       write(dos_file,*) "#    2        Down-spin DOS ", trim(dos_units)
       write(dos_file,*) "#    3        Up-spin Integrated DOS ", trim(intdos_units)
       write(dos_file,*) "#    4        Down-spin Integrated DOS ", trim(intdos_units)
    else
       write(dos_file,*) "#    2        DOS ", trim(dos_units)
       write(dos_file,*) "#    3        Integrated DOS ", trim(intdos_units)
    endif
    write(dos_file, '(1x,a1)') "#"
    write(dos_file, '(1x,a78)') "##############################################################################"

    if(nspins>1) then
       do i=1,nbins
          write(dos_file, *) E(i), dos(i,1), -dos(i,2),  intdos(i,1), -intdos(i,2)
       enddo
    else
       do i=1,nbins
          write(dos_file, *) E(i), dos(i,1), intdos(i,1)
       enddo
    endif
    close(dos_file)
  end subroutine write_dos
  
  
endmodule od_dos
