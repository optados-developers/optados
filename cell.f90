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
 use od_constants,only : dp 
 implicit none

 private
!-------------------------------------------------------------------------!
! R E A L   S P A C E   V A R I A B L E S
 real(kind=dp), public, save :: real_lattice(1:3,1:3)
 real(kind=dp), public, save :: recip_lattice(1:3,1:3)
 real(kind=dp), public, save :: cell_volume
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! R E C I P R O C A L   S P A C E   V A R I A B L E S
 real(kind=dp), allocatable, public, save :: kpoint_r(:,:)
 real(kind=dp), allocatable, public, save :: kpoint_r_cart(:,:)
 real(kind=dp), allocatable, public, save :: kpoint_weight(:)
  integer, allocatable,public, save       :: num_kpoints_on_node(:)

 integer, public, save :: nkpoints 
 integer, public, save :: kpoint_grid_dim(3)
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! G L O B A L L Y   A V A I L A B L E   F U N C T I O N S
 public :: cell_find_MP_grid
 public :: cell_calc_lattice
 public :: cell_report_parameters
!-------------------------------------------------------------------------!

 contains

!=========================================================================!
 subroutine cell_find_MP_grid(kpoints,num_kpts,kpoint_grid_dim)
!=========================================================================!
! Find the MP grid from a set of kpoints                                  !
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
! Written by Andrew Morris from the lindos program             11/10/2010 !
!=========================================================================!
  implicit none

  integer, intent(out) :: kpoint_grid_dim(1:3)
  integer, intent(in)  :: num_kpts
  real(kind=dp),intent(in) :: kpoints(1:3,1:num_kpts)

   call kpoint_density(kpoints(1,:), num_kpts, kpoint_grid_dim(1))
   call kpoint_density(kpoints(2,:), num_kpts, kpoint_grid_dim(2))
   call kpoint_density(kpoints(3,:), num_kpts, kpoint_grid_dim(3))

end subroutine cell_find_MP_grid
!=========================================================================!


!=========================================================================
 subroutine kpoint_density(vector,length,points)
!=========================================================================
! Work out the number of kpoints in a given dimension
! Receives a vector of length, runs though the vector looking for the
! closest distance between points. Then invert this distance to find the
! number of k-points in the given dimension. 
!-------------------------------------------------------------------------
! Arguments: vector (in) : vector of kpoints in a particular dimension                              
!            length (in) : length of vector                  
!            points (out): number of kpoints in the directon of vector
!-------------------------------------------------------------------------
! Parent module variables used: None                                      
!-------------------------------------------------------------------------
! Modules used:  None                                                     
!-------------------------------------------------------------------------
! Key Internal Variables: Described below                                                         
!-------------------------------------------------------------------------
! Necessary conditions: None            
!-------------------------------------------------------------------------
! Known Worries: None                                      !
!-------------------------------------------------------------------------
! Written by Andrew Morris from the LinDOS program             11/10/2010 
!=========================================================================
  implicit none
  integer,       intent(in)  :: length
  real(kind=dp), intent(in)  :: vector(1:length)
  integer,       intent(out) :: points

  real(kind=dp) :: real_points
  real(kind=dp) :: vect(1:length)
  real(kind=dp) :: distance

  logical :: small_found
  integer :: i,j
! THIS IS THE BOMB-PROOF 2nd VERSION
! from LinDOS

! Check that there isn't only one-kpoint
  if(length==1) then
   points=1
   return
  endif

! Check that all of the k-points don't have the same value
  i=1
  do
   if(.not.vector(i)==vector(i+1)) exit
   if((i+1)==length) then ! We haven't exited yet
   ! therefore all must be the same
    points=1
    return
   endif
   i=i+1
 enddo

! Now we have more than one k-point and they're not unique
distance=huge(distance)
do i=1,length
 do j=i+1,length
  if(vector(i)>vector(j)) then
   if((vector(i)-vector(j))<distance) distance=(vector(i)-vector(j))
  elseif(vector(j)>vector(i)) then
    if((vector(j)-vector(i))<distance) distance=(vector(j)-vector(i))
  endif
 enddo
enddo


  real_points=1.0_dp/(distance)
  points=int(real_points+0.5_dp)

end subroutine

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
use od_constants, only : pi,bohr
implicit none

  ! THESE ARE IN BOHR, DON'T GET TRIPPED UP AGAIN!
  real_lattice=real_lattice*bohr

  recip_lattice(1,1)=real_lattice(2,2)*real_lattice(3,3)- &
       real_lattice(3,2)*real_lattice(2,3)
  recip_lattice(2,1)=real_lattice(2,3)*real_lattice(3,1)- &
       real_lattice(3,3)*real_lattice(2,1)
  recip_lattice(3,1)=real_lattice(2,1)*real_lattice(3,2)- &
       real_lattice(3,1)*real_lattice(2,2)
  recip_lattice(1,2)=real_lattice(3,2)*real_lattice(1,3)- &
       real_lattice(1,2)*real_lattice(3,3)
  recip_lattice(2,2)=real_lattice(3,3)*real_lattice(1,1)- &
       real_lattice(1,3)*real_lattice(3,1)
  recip_lattice(3,2)=real_lattice(3,1)*real_lattice(1,2)- &
       real_lattice(1,1)*real_lattice(3,2)
  recip_lattice(1,3)=real_lattice(1,2)*real_lattice(2,3)- &
       real_lattice(2,2)*real_lattice(1,3)
  recip_lattice(2,3)=real_lattice(1,3)*real_lattice(2,1)- &
       real_lattice(2,3)*real_lattice(1,1)
  recip_lattice(3,3)=real_lattice(1,1)*real_lattice(2,2)- &
       real_lattice(2,1)*real_lattice(1,2)

  ! * Calculate cell volume
  cell_volume=real_lattice(1,1)*recip_lattice(1,1) + &
       real_lattice(2,1)*recip_lattice(1,2) + &
       real_lattice(3,1)*recip_lattice(1,3)


  if(cell_volume<0.0_dp) then ! Left handed set
    cell_volume=-cell_volume
  endif

  ! Scale reciprocal lattice by 2*pi/volume
  recip_lattice(:,:)=recip_lattice(:,:)*pi*2.0_dp/cell_volume

 end subroutine
 
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
use od_io, only : stdout

  implicit none

  integer :: i

  write(stdout,'(30x,a21)') 'Lattice Vectors (Ang)' 
  write(stdout,101) 'a_1',(real_lattice(1,I), i=1,3)
  write(stdout,101) 'a_2',(real_lattice(2,I), i=1,3)
  write(stdout,101) 'a_3',(real_lattice(3,I), i=1,3)

  write(stdout,*)   
   
  write(stdout,'(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
  write(stdout,101) 'b_1',(recip_lattice(1,I), i=1,3)
  write(stdout,101) 'b_2',(recip_lattice(2,I), i=1,3)
  write(stdout,101) 'b_3',(recip_lattice(3,I), i=1,3)

  write(stdout,*) 
  write(stdout,'(19x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume
 
  write(stdout,'(2x,a7)') '(Ang^3)'
  write(stdout,*) 
  
  return
101 format(20x,a3,2x,3F11.6)
  
 end subroutine cell_report_parameters

endmodule od_cell
