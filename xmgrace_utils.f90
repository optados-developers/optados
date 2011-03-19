!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
module xmgrace_utils
 ! use od_constants, only : dp

  implicit none

  public :: xmgu_setup
  public :: xmgu_legend
  public :: xmgu_title
  public :: xmgu_subtitle
  public :: xmgu_axis
  public :: xmgu_data
  
  private



contains

  !==================================================================!
  subroutine xmgu_setup(unit)
  !==================================================================!
    implicit none
    integer, intent(in)      :: unit
   
     write(unit,*) 'autoscale'
     write(unit,*) 'world xmin 0.0'
     write(unit,*) 'world xmax 1.0'
     write(unit,*) 'world ymax 0.1' 
     
     write(unit,*) 'frame linestyle 1'
     write(unit,*) 'frame linewidth 2.5' 
     write(unit,*) 'frame color 1'
     write(unit,*) 'frame pattern 1'
     write(unit,*) 'frame background color 0'
     write(unit,*) 'frame background pattern 0'
     
     write(unit,*) 'altxaxis  off'
     write(unit,*) 'altyaxis  off'


  end subroutine xmgu_setup
  !==================================================================!


 !==================================================================!
  subroutine xmgu_legend(unit)
   !==================================================================!
    implicit none
    integer, intent(in)      :: unit
    
    ! Clearly this will not be true
     write(unit,*) 'legend off'
    
  end subroutine xmgu_legend  
   !==================================================================!
   

  !==================================================================!
  subroutine xmgu_title(unit,title)
  !==================================================================!
    implicit none
    
    integer, intent(in)      :: unit
    character(80),intent(in) :: title
    
    write(unit,*) 'title "'//trim(title)//'"'
    write(unit,*) "title font 4"
    write(unit,*) "title size 1.500000"
    write(unit,*) "title color 1"
   
  end subroutine xmgu_title
  !==================================================================!
  
  
  !==================================================================!
  subroutine xmgu_subtitle(unit,subtitle)
  !==================================================================!
    implicit none
    
    integer, intent(in)      :: unit
    character(80),intent(in) :: subtitle
    
    write(unit,*) 'subtitle '//trim(subtitle)//'"'
    write(unit,*) "subtitle font 4"
    write(unit,*) "subtitle size 1.000000"
    write(unit,*) "subtitle color 1"

  
  end subroutine xmgu_subtitle
  !==================================================================!
  
  
  !==================================================================!
  subroutine xmgu_axis(unit,axis,label)
  !==================================================================!
    implicit none
    
    integer, intent(in)      :: unit
    character(80),intent(in) :: label
    character(1),intent(in)  :: axis
    
    character(5) :: axis_name
    
    if(axis=="x") then
     axis_name="xaxis"
    elseif(axis=="y") then
     axis_name="yaxis"
    endif
    
    write(unit,*) trim(axis_name)//" type zero true"
    write(unit,*) trim(axis_name)//" bar on"   
    write(unit,*) trim(axis_name)//" bar color 1"
    write(unit,*) trim(axis_name)//" bar linestyle 1"
    write(unit,*) trim(axis_name)//" bar linewidth 2.5"
    write(unit,*) trim(axis_name)//' label "'//trim(label)//'"'
    write(unit,*) trim(axis_name)//" label layout para"
    write(unit,*) trim(axis_name)//" label place auto" 
    write(unit,*) trim(axis_name)//" label char size 1.000000"
    write(unit,*) trim(axis_name)//" label font 4"
    write(unit,*) trim(axis_name)//" label color 1"
    write(unit,*) trim(axis_name)//" label place normal" 
    write(unit,*) trim(axis_name)//" tick on"
    write(unit,*) trim(axis_name)//" tick major 0.2"
    write(unit,*) trim(axis_name)//" tick minor ticks 1"
    write(unit,*) trim(axis_name)//" tick default 6"
    write(unit,*) trim(axis_name)//" tick place rounded true"
    write(unit,*) trim(axis_name)//" tick in"
    write(unit,*) trim(axis_name)//" tick major size 1.000000"
    write(unit,*) trim(axis_name)//" tick major color 1"
    write(unit,*) trim(axis_name)//" tick major linewidth 2.5"
    write(unit,*) trim(axis_name)//" tick major linestyle 1"
    write(unit,*) trim(axis_name)//" tick major grid off"
    write(unit,*) trim(axis_name)//" tick minor color 1"
    write(unit,*) trim(axis_name)//" tick minor linewidth 2.5"
    write(unit,*) trim(axis_name)//" tick minor linestyle 1"
    write(unit,*) trim(axis_name)//" tick minor grid off"
    write(unit,*) trim(axis_name)//" tick minor size 0.500000"
    write(unit,*) trim(axis_name)//" ticklabel on"
    write(unit,*) trim(axis_name)//" ticklabel format general"
    write(unit,*) trim(axis_name)//" ticklabel prec 5"
    write(unit,*) trim(axis_name)//" ticklabel start type spec"
    write(unit,*) trim(axis_name)//" ticklabel start 0.200000"
    write(unit,*) trim(axis_name)//" ticklabel stop type spec"
    write(unit,*) trim(axis_name)//" ticklabel stop 0.800000"
    write(unit,*) trim(axis_name)//" ticklabel font 4"

  end subroutine xmgu_axis
  !==================================================================!
  
  
    
  !==================================================================!
  subroutine xmgu_data(unit)
  !==================================================================!
    integer, intent(in)      :: unit
 ! s0 hidden false
!s0 type xy
!s0 symbol 9
!s0 symbol size 1.000000
!s0 symbol color 4
!s0 symbol pattern 1
!s0 symbol fill color 4
!s0 symbol fill pattern 1
!s0 symbol linewidth 3.0
!s0 symbol linestyle 1
!s0 symbol char 65
!s0 symbol char font 0
!s0 symbol skip 0
!s0 line type 0
!s0 line linestyle 1
!s0 line linewidth 1.0
!s0 legend  ""
!s0 avalue on
!s0 avalue type 4
!s0 avalue char size 0.80000
!s0 avalue font 4
!s0 avalue color 1
!s0 avalue rot 0
!s0 avalue format general
!s0 avalue prec 4
!s0 avalue prepend ""
!s0 avalue append ""
!s0 avalue offset 0.000000 , 0.01000
  end subroutine xmgu_data



  	
end module xmgrace_utils