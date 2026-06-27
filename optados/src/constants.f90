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

  character(len=6), parameter, public :: optados_version = " 1.3 "
  character(len=14), parameter, public :: copyright = " (c) 2010-2025"

  integer, parameter, public :: dp = selected_real_kind(15, 300)

  real(kind=dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp
  real(kind=dp), parameter, public :: H2eV = 27.21138342902473_dp
  real(kind=dp), parameter, public :: bohr2ang = 0.52917720859_dp
  real(kind=dp), parameter, public :: inv_sqrt_two_pi = 0.3989422804014326779399460599_dp
  real(kind=dp), parameter, public :: twopi = 6.283185307179586476925286766559005768394_dp
  real(kind=dp), parameter, public :: sqrt_two = 1.414213562373095048801688724209698079_dp
  complex(dp), parameter, public   :: cmplx_0 = (0.0_dp, 0.0_dp)
  complex(dp), parameter, public   :: cmplx_i = (0.0_dp, 1.0_dp)

  ! Used in phonon_eels.
  real(kind=dp), parameter, public :: amu_to_me = 1822.888486209_dp
  real(kind=dp), parameter, public :: speed_of_light_SI = 8.8541878188E-12
  real(kind=dp), parameter, public :: epsilon_0_SI = 8.8541878188E-12
  real(kind=dp), parameter, public :: elec_charge_SI = 1.602176634E-19_dp
  real(kind=dp), parameter, public :: hbar_j_s = 1.054571817e-34_dp
  real(kind=dp), parameter, public :: electron_mass_kg = 9.1093837015e-31_dp
  real(kind=dp), parameter, public :: eV_to_joule = 1.602176634e-19_dp
  real(kind=dp), parameter, public :: eV_to_hartree = 0.03674932217565499_dp
  real(kind=dp), parameter, public :: meV_to_hartree = 3.674932217565499E-5_dp
  real(kind=dp), parameter, public :: inv_cm_to_meV = 0.12398419843320026_dp
  real(kind=dp), parameter, public :: ang2bohr = 1.8897259886_dp
  real(kind=dp), parameter, public :: nanometre2bohr = 18.89726124565062_dp
  real(kind=dp), parameter, public :: fourpi = 12.56637061435917295385057353311801_dp
  real(kind=dp), parameter, public :: fine_structure_constant = 0.0072973525643_dp
  real(kind=dp), parameter, public :: meV_per_K = 0.08617333262145_dp

  character(len=3), parameter, dimension(109), public :: periodic_table_name = (/ &
       & 'H ', 'He', &
       & 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
       & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
       & 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
       & 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
       & 'Cs', 'Ba', &
       & 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
       & 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
       & 'Fr', 'Ra', &
       & 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', &
       & 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt'/)

  real(kind=dp), parameter, dimension(109), public :: preset_debye_waller = 0.0_dp

end module od_constants
