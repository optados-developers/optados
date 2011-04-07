!=========================================================================!
! Module: Constants                                                       !
! For global constants                                                    !
!-------------------------------------------------------------------------!
! Modules used:  None                                                     !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
!-------------------------------------------------------------------------!
! Written by Andrew Morris                                     11/10/2010 !
!=========================================================================!
module od_constants
  implicit none
  
  private ! unless otherise stated

 character(len=6), parameter, public :: optados_version=" 0.9  "

 integer, parameter, public :: dp=selected_real_kind(15,300)

 real(kind=dp), parameter, public :: pi=3.141592653589793238462643383279502884197_dp
 real(kind=dp), parameter, public :: H2eV=27.21138342902473_dp
 real(kind=dp), parameter, public :: bohr2ang=0.52917720859_dp
 real(kind=dp), parameter, public :: inv_sqrt_two_pi=0.3989422804014326779399460599_dp
 real(kind=dp), parameter, public :: twopi=6.283185307179586476925286766559005768394_dp
 real(kind=dp), parameter, public :: sqrt_two=1.414213562373095048801688724209698079_dp
 complex(dp),   parameter, public :: cmplx_0=(0.0_dp, 0.0_dp)
 complex(dp),   parameter, public :: cmplx_i=(0.0_dp, 1.0_dp)

endmodule od_constants
