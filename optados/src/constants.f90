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
  character(len=14), parameter, public :: copyright = " (c) 2010-2022"

  integer, parameter, public :: dp = selected_real_kind(15, 300)

  real(kind=dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp
  real(kind=dp), parameter, public :: H2eV = 27.21138342902473_dp
  real(kind=dp), parameter, public :: bohr2ang = 0.52917720859_dp
  real(kind=dp), parameter, public :: inv_sqrt_two_pi = 0.3989422804014326779399460599_dp
  real(kind=dp), parameter, public :: twopi = 6.283185307179586476925286766559005768394_dp
  real(kind=dp), parameter, public :: sqrt_two = 1.414213562373095048801688724209698079_dp
  complex(dp), parameter, public :: cmplx_0 = (0.0_dp, 0.0_dp)
  complex(dp), parameter, public :: cmplx_i = (0.0_dp, 1.0_dp)

  !Optics constants
  real(kind=dp), parameter, public :: epsilon_0 = 8.8541878176E-12_dp
  real(kind=dp), parameter, public :: e_charge = 1.602176487E-19_dp
  real(kind=dp), parameter, public :: e_mass = 9.10938215E-31_dp
  real(kind=dp), parameter, public :: hbar = 1.054571628E-34_dp
  real(kind=dp), parameter, public :: c_speed = 299792458.0_dp

  !Photoemission constants
  real(kind=dp), parameter, public :: epsilon_zero = 55.26349406_dp     !e^2 GeV^-1 fm^-1
  real(kind=dp), parameter, public :: j_to_ev = 6.24150934E+18_dp      !J eV^-1
  real(kind=dp), parameter, public :: ev_to_j = 1.602176565E-19_dp          !eV J^-1
  real(kind=dp), parameter, public :: ev_to_hartree = 0.03674932379_dp          !eV Ha^-1
  real(kind=dp), parameter, public :: rad_to_deg = 57.2957795_dp          !deg rad^-1
  real(kind=dp), parameter, public :: deg_to_rad = 0.0174532925_dp      !rad^-1 to deg
  real(kind=dp), parameter, public :: boltzmann = 1.38064852E-23_dp      !J K^-1
  real(kind=dp), parameter, public :: kB = 8.617333262E-5_dp                ! ev K^-1
  ! Constants for field emission in Photoemission Module
  real(kind=dp), parameter, public :: b_factor = 0.6830890_dp            ! eV^(-3/2) V A^-1
  real(kind=dp), parameter, public :: p1 = 0.03270530446_dp
  real(kind=dp), parameter, public :: p2 = 0.009157798739_dp
  real(kind=dp), parameter, public :: p3 = 0.002644272807_dp
  real(kind=dp), parameter, public :: p4 = 0.00008987173811_dp
  real(kind=dp), parameter, public :: q1 = 0.1874993441_dp
  real(kind=dp), parameter, public :: q2 = 0.01750636947_dp
  real(kind=dp), parameter, public :: q3 = 0.005527069444_dp
  real(kind=dp), parameter, public :: q4 = 0.001023904180_dp

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

  ! Based on: S. Alvarez, Dalton Transactions, 2013, 42, 8617-8636.
  ! If value == 1.0_dp then no value was supplied, subject to change
  real(kind=dp), parameter, dimension(109), public :: periodic_table_vdw = (/ &
      & 1.20_dp, 1.43_dp, &
      & 2.12_dp, 1.98_dp, 1.91_dp, 1.77_dp, 1.66_dp, 1.50_dp, 1.46_dp, 1.58_dp, &
      & 2.50_dp, 2.51_dp, 2.25_dp, 2.19_dp, 1.90_dp, 1.89_dp, 1.82_dp, 1.83_dp, &
      & 2.73_dp, 2.62_dp, 2.58_dp, 2.46_dp, 2.42_dp, 2.45_dp, 2.45_dp, 2.44_dp, 2.40_dp, 2.40_dp, 2.38_dp, 2.39_dp, 2.32_dp, &
      & 2.29_dp, 1.88_dp, 1.82_dp, 1.86_dp, 2.25_dp, &
      & 3.21_dp, 2.84_dp, 2.75_dp, 2.52_dp, 2.56_dp, 2.45_dp, 2.44_dp, 2.46_dp, 2.44_dp, 2.15_dp, 2.53_dp, 2.49_dp, 2.43_dp, &
      & 2.42_dp, 2.47_dp, 1.99_dp, 2.04_dp, 2.06_dp, &
      & 3.48_dp, 3.03_dp, &
      & 2.98_dp, 2.88_dp, 2.92_dp, 2.95_dp, 1.0_dp, 2.90_dp, 2.87_dp, 2.83_dp, 2.79_dp, 2.87_dp, 2.81_dp, 2.83_dp, 2.79_dp, &
      & 2.80_dp, 2.74_dp, &
      & 2.63_dp, 2.53_dp, 2.57_dp, 2.49_dp, 2.48_dp, 2.41_dp, 2.29_dp, 2.32_dp, 2.45_dp, 2.47_dp, 2.60_dp, 2.54_dp, 1.0_dp, &
      & 1.0_dp, 1.0_dp, &
      & 1.0_dp, 1.0_dp, &
      & 2.8_dp, 2.93_dp, 2.88_dp, 2.71_dp, 2.82_dp, 2.81_dp, 2.83_dp, 3.05_dp, 3.4_dp, 3.05_dp, 2.7_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
      & 1.0_dp, &
      & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/)

end module od_constants
