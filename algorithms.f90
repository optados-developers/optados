!=========================================================================!
! Module: algorithms                                                      !
! For low-level algorithms that are not specific to electronic structure  !
!-------------------------------------------------------------------------!
! Modules used:  constants - dp. inv_sqrt_two_pi                          !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below                                                         !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris and Chris Picakrd (so far)          11/10/2010 !
!=========================================================================!
module od_algorithms
  use od_constants, only : dp, inv_sqrt_two_pi
  implicit none
 
  private

  public :: gaussian
  public :: heap_sort
  public :: utility_lowercase
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart


 contains

!=========================================================================!
  function gaussian(m,w,x)
!=========================================================================!
! ** Return value of Gaussian(mean=m,width=w) at position x
! I don't know who's this function originally was, CJP? MIJP?
!=========================================================================!
    implicit none

    real(kind=dp), intent(in) :: m,w,x
    real(kind=dp)             :: gaussian

    if(0.5_dp*((x-m)/w)**2.gt.30.0_dp) then

    gaussian=0.0_dp
    return
    else
    gaussian=inv_sqrt_two_pi*exp(-0.5_dp*((x-m)/w)**2)/w
    endif
    return
  end function gaussian


!=========================================================================! 
  subroutine heap_sort(num_items,weight)
 !=========================================================================!
 ! This subroutine sorts the list of weights into descending order.        !
 ! On exit, if present, the array of indexes contains the original index   !
 ! of each item.                                                           !
 !                                                                         !
 ! This is a heap sort                                                     !
 !-------------------------------------------------------------------------!
 ! Arguments:                                                              !
 !   num_items (input) :: The number of items to sort                      !
 !   weight (in/out) :: The weights of each item. On exit these are        !
 !                      sorted into descending order.                      !
 !-------------------------------------------------------------------------!
 ! Parent module variables used:                                           !
 !   None                                                                  !
 !-------------------------------------------------------------------------!
 ! Modules used:                                                           !
 !   None                                                                  !
 !-------------------------------------------------------------------------!
 ! Key Internal Variables:                                                 !
 !   None                                                                  !
 !-------------------------------------------------------------------------!
 ! Necessary conditions:                                                   !
 !   None                                                                  !
 !-------------------------------------------------------------------------!
 ! Written by Chris Pickard 22nd May 2009                                  !
 !=========================================================================!
   implicit none

   ! Arguments

   integer, intent(in) :: num_items
   real(kind=dp), dimension(num_items), intent(inout) :: weight

   ! Local variables

   integer :: i,ir,j,l ! Loop counters
   real(kind=dp) :: wta

   if(num_items.lt.2) return

   l=num_items/2+1
   ir=num_items

   do
      if(l.gt.1) then
         l=l-1
         wta=weight(l)
      else
         wta=weight(ir)
         weight(ir)=weight(1)
         ir=ir-1
         if(ir.eq.1) then
            weight(1)=wta
            return
         end if
      end if
      i=l
      j=l+l
20     if(j.le.ir) then
         if(j.lt.ir) then
            if(weight(j).lt.weight(j+1)) j=j+1
         end if
         if(wta.lt.weight(j)) then
            weight(i)=weight(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         goto 20
      end if
      weight(i)=wta
   end do

 end subroutine heap_sort

!=========================================================================!
 function utility_lowercase(string)
!=========================================================================!
! Takes a string and converts to lowercase characters       
!=========================================================================!
     
    use od_io, only : maxlen
    implicit none
  
    character(len=*), intent(in) :: string
    character(len=maxlen)        :: utility_lowercase

    integer :: iA,iZ,idiff,ipos,ilett

    iA = ichar('A')
    iZ = ichar('Z')
    idiff = iZ-ichar('z')

    utility_lowercase = string

    do ipos=1,len(string)
       ilett = ichar(string(ipos:ipos))
       if ((ilett.ge.iA).and.(ilett.le.iZ)) &
            utility_lowercase(ipos:ipos)=char(ilett-idiff)
    enddo

    utility_lowercase = trim(adjustl(utility_lowercase))

    return

  end function utility_lowercase


    !===================================================================
         subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(1,i)*frac(1) + real_lat(2,i)*frac(2) + real_lat(3,i)*frac(3) 
    end do

    return

  end subroutine utility_frac_to_cart

    !===================================================================
         subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    use od_constants, only : twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3) 
    end do

    frac=frac/twopi


    return

  end subroutine utility_cart_to_frac



endmodule od_algorithms
