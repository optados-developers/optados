module dos
 use constants, only : dp
 implicit none

! type bands_info
!   real(kind=dp) :: energies
!   real(kind=dp) :: gradients   
!!   real(kind=dp) :: curvatures
! end type bands_info
!
! type kpoints_info
!   real(kind=dp) :: r(1:3)
!   real(kind=dp) :: weights 
! end type kpoints_info
!
!type dos_info
!   real(kind=dp) :: gaussion
!   real(kind=dp) :: linear
!   real(kind=dp) :: smear  
! end type kpoints_info

! private ! unless otherwise indicated

!type(bands_info)  , allocatable, public :: band(:)
!type(dos_info)    , allocatable, public :: int_dos(:)
!type(dos_info)    , allocatable, public :: dos(:), weighted_dos(:)
!type(kpoints_info), allocatable, public :: kpoint(:)



!contains
! subroutine dos_slice
! end subroutine dos_slice
!
! subroutine dos_calculate(band,dos,weights,weighted_dos)
!  implicit none
!  type(bands_info)  , allocatable, public :: band(:), band_slice(:)
!  type(dos_info)    , allocatable, public :: int_dos(:), int_dos_slice(:)
!  type(dos_info)    , allocatable, public :: dos(:), dos_slice(:), weighted_dos(:), weighted_slice(:)
!  type(kpoints_info), allocatable, public :: kpoint(:), kpoint_slice(:)
!
!  call dos_calculate_slice
! 
!  call dos_merge
!
!eend subroutine dos_calculate
!
! subroutine dos_calcualte_slice
! end subroutine dos_calculate_slice
!
! subroutine dos_merge
! end subroutine dos_merge

 ! GOOD ISN'T IT?

endmodule dos
