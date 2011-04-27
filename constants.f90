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


  character(len=3), parameter, dimension(109), public :: periodic_table_name= (/ &
       & 'H ',                                                                                'He', &
       & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
       & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
       & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
       & 'Cs','Ba', &
       & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
       & 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
       & 'Fr','Ra', &
       & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr', &
       & 'Rf','Db','Sg','Bh','Hs','Mt' /)



endmodule od_constants
