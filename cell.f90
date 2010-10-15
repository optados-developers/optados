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
module cell
 use constants, only : dp
 implicit none

 private

 public find_MP_grid

 contains

 subroutine find_MP_grid(kpoints,num_kpts,nx,ny,nz)
!=========================================================================!
! Find the MP grid from a set of kpoints                                  !
!-------------------------------------------------------------------------!
! Arguments: kpoints - an array of kpoints                                !
!            num_kpts - size of the kpoint array                          !
!-------------------------------------------------------------------------!
! Returns: nx,ny,nz - the number of kpoints in each dimension             !
!-------------------------------------------------------------------------!
! Parent module variables used: None                                      !
!-------------------------------------------------------------------------!
! Modules used:  None                                                     !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below                                                         !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris from the lindos program             11/10/2010 !
!=========================================================================!
  implicit none

  integer, intent(out) :: nx,ny,nz
  integer, intent(in)  :: num_kpts
  real(kind=dp),intent(in) :: kpoints(1:3,1:num_kpts)

   call kpoint_density(kpoints(1,:), num_kpts, nx)
   call kpoint_density(kpoints(2,:), num_kpts, ny)
   call kpoint_density(kpoints(3,:), num_kpts, nz)
 contains

 subroutine kpoint_density(vector,length,points)
  implicit none
  integer, intent(in)        :: length
  real(kind=dp), intent(in)  :: vector(1:length)
  real(kind=dp)   :: real_points
  real(kind=dp)   :: vect(1:length)
  integer, intent(out) :: points
  real(kind=dp) :: distance

  logical :: small_found
  integer :: i,j
! THIS IS THE BOMB-PROOF 2nd VERSION


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

end subroutine

endmodule cell
