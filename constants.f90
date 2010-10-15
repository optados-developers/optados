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
module constants
  implicit none
  
  private ! unless otherise stated

 integer, parameter, public :: dp=selected_real_kind(15,300)

 real(kind=dp), parameter, public :: pi=3.141592653589793238462643383279502884197_dp
 real(kind=dp), parameter, public :: H2eV=27.21138342902473_dp
 real(kind=dp), parameter, public :: bohr=0.52917720859_dp
 real(kind=dp), parameter, public :: inv_sqrt_two_pi=0.3989422804014326779399460599_dp


endmodule constants
