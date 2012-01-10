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
module od_dos
  use od_constants, only : dp
  !use od_dos_utils, only : dos_utils_calculate
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
    use od_io,        only : stdout,io_time
    use od_dos_utils, only : E, dos_fixed, intdos_fixed, dos_adaptive, &
    &intdos_adaptive, dos_linear, intdos_linear,dos_utils_calculate
    use od_parameters,only : fixed,adaptive,linear
    use od_comms,     only : on_root
    use od_electronic,only : nspins,efermi
       
       
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
    use od_parameters, only : dos_nbins, dos_per_volume, output_format
    use od_io,         only : seedname, io_file_unit,io_date,io_error, stdout

    implicit none
    real(dp), intent(in) :: E(dos_nbins)
    real(dp), intent(in) :: dos(dos_nbins,nspins)
    real(dp), intent(in) :: intdos(dos_nbins,nspins)
    character(len=*), intent(in) :: dos_name
    integer :: i, dos_file, ierr
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=22) :: dos_units, intdos_units


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
       do i=1,dos_nbins
          write(dos_file, *) E(i), dos(i,1), -dos(i,2),  intdos(i,1), -intdos(i,2)
       enddo
    else
       do i=1,dos_nbins
          write(dos_file, *) E(i), dos(i,1), intdos(i,1)
       enddo
    endif
    close(dos_file)

   if(trim(output_format)=="xmgrace") then
     call write_dos_xmgrace(dos_name,E,dos)
   elseif(trim(output_format)=="gnuplot") then 
     write(stdout,*)  " WARNING: GNUPLOT output not yet available, calling xmgrace"
     call write_dos_xmgrace(dos_name,E,dos)
!     call write_dos_gnuplot(dos_name,E,dos)
   else
     write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
  endif



  end subroutine write_dos

    !=============================================================================== 
  subroutine write_dos_xmgrace(dos_name,E,dos)
    !=============================================================================== 
   use xmgrace_utils
   use od_parameters, only : dos_nbins
   use od_electronic, only : nspins
   use od_io,         only : io_file_unit,io_error,seedname 
   implicit none 

   real(dp), intent(in) :: E(dos_nbins)
    real(dp), intent(in) :: dos(dos_nbins,nspins)
  
   real(dp) :: min_x, max_x, min_y, max_y
  
   integer :: batch_file,ierr
   character(len=*), intent(in) :: dos_name

   batch_file=io_file_unit()
   open(unit=batch_file,file=trim(seedname)//'.'//trim(dos_name)//'.agr',iostat=ierr)
   if(ierr.ne.0) call io_error(" ERROR: Cannot open xmgrace batch file in dos: write_dos_xmgrace")

   min_x=minval(E)
   max_x=maxval(E)
   
   min_y=0
   max_y=maxval(dos)
   if(nspins>1) then
     min_y=-max_y
   endif


   call  xmgu_setup(batch_file)
   call  xmgu_legend(batch_file)
   call  xmgu_title(batch_file, min_x, max_x, min_y, max_y, "Electronic Density of States")
   call  xmgu_subtitle(batch_file,"Generated by OptaDOS")
    
   call  xmgu_axis(batch_file,"x","Energy eV")
   call  xmgu_axis(batch_file,"y","eDOS")

   if(nspins>1) then
      call xmgu_data_header(batch_file,0,1,"up-spin channel")
      call xmgu_data_header(batch_file,1,2,"down-spin channel")
      call xmgu_data(batch_file,0,E(:),dos(:,1))
      call xmgu_data(batch_file,1,E(:),-dos(:,2))
   else
      call xmgu_data_header(batch_file,0,1,"Total DOS")
      call xmgu_data(batch_file,0,E(:),dos(:,1))  
   endif
  
  
  close(batch_file)

  end subroutine write_dos_xmgrace
  
  
endmodule od_dos
