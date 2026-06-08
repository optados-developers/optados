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
!===========================================================================!
! MODULE od_phonon_eels                                                     !
!                                                                           !
! J. A. J. Whaley-Baldwin, April 2026                                       !
!                                                                           !
! This module contains routines for calculating a phonon EELS spectrum,     !
! given phonon data from a CASTEP .phonon file.                             !
!                                                                           !
! Both impact and aloof methods are implemented.                            !
!===========================================================================!
module od_phonon_eels

  ! OptaDOS comms.
  use od_comms, only: on_root

  ! OptaDOS IO.
  use od_io, only: seedname, stdout, io_error, maxlen, io_file_unit

  ! OptaDOS constants.
  use od_constants, only: dp, pi, twopi, fourpi, cmplx_i, cmplx_0, amu_to_me, ang2bohr, nanometre2bohr, &
  fine_structure_constant, inv_cm_to_meV, eV_to_hartree, meV_to_hartree, meV_per_K

  ! OptaDOS parameters.
  use od_parameters, only: phonon_eels_task, vibeels_min_energy, vibeels_max_energy, vibeels_spacing, &
  vibeels_impact_broadening, vibeels_use_electronic_affs, vibeels_aloof_intrinsic_broadening, vibeels_aloof_loss_broadening, &
  vibeels_aloof_surface_plane, vibeels_aloof_beam_direction, vibeels_aloof_include_offdiag, vibeels_aloof_electron_beam_energy, &
  vibeels_aloof_impact_parameter, vibeels_aloof_phi_spacing, vibeels_reorder_phonon_bands, iprint

  ! Lookup tables for atomic quantities.
  use od_phonon_eels_tables, only: get_atomic_number, get_cromer_mann_coeff, get_peng_electron_coeff

  ! Essential data for vib-EELS.
  use od_read_odd, only: parsed_eels_data_from_odd, odd_read_vib_eels_data, inf_dielectric_tensor, born_eff_ch_tensor, &
  partial_charges, atomic_displacement_params, ADP_temperature

  implicit none

  private

  ! Exposed routines.
  public :: phonon_eels_calculate

  ! ******************************************  vib-EELS QUANTITIES  ********************************************** !

  ! These arrays will store the calculated vib-EELS quantities.
  ! They are allocated and populated by calling the various subroutines in this module.
  complex(kind=dp),allocatable,save  :: mode_resolved_born_eff_ch_tensor(:,:)          ! (mode_idx,i)
  complex(kind=dp),allocatable,save  :: osc_strength_tensor(:,:,:)                     ! (mode,i,j)
  complex(kind=dp),allocatable,save  :: lf_dielectric_tensor(:,:,:)                    ! (i,j,omega)
  complex(kind=dp),allocatable,save  :: polarizability(:,:)                            ! (omega,phi)
  complex(kind=dp),allocatable,save  :: aloof_loss_probability(:)                      ! (omega)
  real(kind=dp),allocatable,save     :: phonon_occupations(:,:)                        ! (mode_idx,qpt_idx)
  real(kind=dp),allocatable,save     :: eels_intensity(:,:)                            ! (mode_idx,qpt_idx)
  real(kind=dp),allocatable          :: impact_heatmap(:,:)                            ! (qpt_idx,omega)

  ! *******************************************  DATA TO BE PARSED IN  ******************************************** !

  ! Basic data that is currently parsed in from the .phonon file.
  ! These could also be parsed in from a .cell file.
  integer,save,public                       :: n_ions
  real(kind=dp),save,public                 :: real_lattice(3,3)                       ! (v,i)
  real(kind=dp),allocatable,save,public     :: atomic_positions(:,:)                   ! (i,atom_idx)
  character(len=2),allocatable,save,public  :: atomic_species(:)                       ! (atom_idx)

  ! These are parsed in uniquely from the .phonon file.
  integer,save,public                       :: n_branches
  integer,save,public                       :: n_qpts
  real(kind=dp),allocatable,save,public     :: atomic_masses(:)                        ! (atom_idx)
  real(kind=dp),allocatable,save,public     :: qpoint_positions(:,:)                   ! (i,qpt_idx)
  real(kind=dp),allocatable,public          :: qpoint_weights(:)                       ! (qpt_idx)
  complex(kind=dp),allocatable,save,public  :: phonon_eigenvectors(:,:,:,:)            ! (mode_idx,atom_idx,dir,qpt_idx)
  real(kind=dp),allocatable,save,public     :: phonon_eigenvalues(:,:)                 ! (mode_idx,qpt_idx)
  complex(kind=dp),allocatable,save,public  :: gamma_eigenvectors(:,:,:)               ! (mode_idx,atom_idx,dir)
  real(kind=dp),allocatable,save,public     :: gamma_eigenvalues(:)                    ! (mode_idx)

  ! **********************************************  OTHER DATA  **************************************************** !

  ! These are set after all the data has been parsed in.
  integer,save,public                       :: N_energies
  integer,save,public                       :: N_phi_values
  real(kind=dp),allocatable,save,public     :: energies(:)                             ! Energies (in meV) to evaluate EELS quantities
  real(kind=dp),allocatable,save,public     :: omegas_hartree(:)                       ! Omegas (in Hartree) to evaluate EELS quantities
  real(kind=dp),allocatable,save,public     :: phis(:)                                 ! Electron EELS phi values (radian)
  real(kind=dp),save,public                 :: recip_lattice(3,3)                      ! (v,i)
  real(kind=dp),save,public                 :: cell_volume                             ! Volume of the real-space unit cell (Ang^3)
  integer,save,public                       :: N_hsps                                  ! Number of high-symmetry points along path
  integer,allocatable,save,public           :: hsp_idxs(:)                             ! Indexes of the high-symmetry points along path
  real(kind=dp),allocatable,save,public     :: path_spacing_norms(:)                   ! q-point spacing along path (inverse Angstrom)

  ! *************************************************  FLAGS  ****************************************************** !

  ! Flag to check whether phonon data has been successfully parsed.
  logical,save,public                       :: parsed_phonon_file      = .false.

  ! Flag to check whether high-symmetry points along q-point path have been detected.
  logical,save,public                       :: detected_hsps           = .false.

  ! Flag to notify whether a rotation was applied to orient the crystal surface for the aloof calculations, or not.
  logical,save,public                       :: aloof_rotation_applied  = .false.

  ! Flag to check whether all relevant prerequisite data has been parsed, and that frequency/phi arrays have been set.
  logical,save,public                       :: phonon_eels_prep_done   = .false.

  ! Flags to check whether various quantities have been calculated, or not.
  logical,save                              :: calculated_mode_resolved_bec_tensor = .false.
  logical,save                              :: calculated_osc_strength_tensor      = .false.
  logical,save                              :: calculated_lf_eps                   = .false.
  logical,save                              :: calculated_polarizability           = .false.
  logical,save                              :: calculated_aloof_loss_probability  = .false.
  logical,save                              :: calculated_thermal_occupations      = .false.
  logical,save                              :: calculated_eels_intensity           = .false.
  logical,save                              :: calculated_impact_heatmap           = .false.

  ! ************************************************  CONTROL  ****************************************************** !

  ! Output file writing.
  logical :: write_eps_lf_to_file                  = .true.
  logical :: write_alpha_lf_to_file                = .true.
  logical :: write_osc_strength_tensor_to_file     = .true.
  logical :: write_eels_intensity_to_file          = .true.
  logical :: write_impact_heatmap_to_file          = .true.
  logical :: write_aloof_loss_probability_to_file = .true.

  ! Enable / Disable optional features.
  logical :: eps_lf_skip_acoustic_modes = .true.        ! Whether to explicitly skip the acoustic modes in the calculation of eps_lf
  logical :: alpha_enforce_analytic_continuity = .true. ! Whether to force analytic continuity in the calculation of alpha

  !******************************************************************************************************************!

  contains

    !===========================================================================!
    subroutine phonon_eels_calculate
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                         !
    !                                                                           !
    ! This subroutine wraps all of the functionality contained within the       !
    ! phonon_eels module. Depending on the value of phonon_eels_task, it calls  !
    ! either impact or aloof subroutines, or both if this is set to 'all'.      !
    !===========================================================================!
      implicit none

      ! No parallelism in phonon_eels at present...
      if (on_root) then

        if ( index(phonon_eels_task, 'impact') > 0 ) then
          CALL phonon_eels_impact_calculate_intensity
          CALL phonon_eels_impact_calculate_heatmap
          CALL phonon_eels_write_affs_and_dwfs_to_file

        else if ( index(phonon_eels_task, 'aloof') > 0 ) then
          CALL phonon_eels_calculate_lf_dielectric_tensor
          CALL phonon_eels_aloof_calculate_polarizability
          CALL phonon_eels_aloof_calculate_loss_probability

        else if ( index(phonon_eels_task, 'all') > 0 ) then
          CALL phonon_eels_impact_calculate_intensity
          CALL phonon_eels_impact_calculate_heatmap
          CALL phonon_eels_write_affs_and_dwfs_to_file
          CALL phonon_eels_calculate_lf_dielectric_tensor
          CALL phonon_eels_aloof_calculate_polarizability
          CALL phonon_eels_aloof_calculate_loss_probability

        else
          CALL io_error("ERROR: phonon_eels_task not recognized. Should be one of 'impact', 'aloof', or 'all'")

        end if

      end if

    end subroutine phonon_eels_calculate

    !===========================================================================!
    !                          *** IMPACT ROUTINES ***                          !
    !                                                                           !
    ! The routines here are used for the impact method, as described in         !
    ! Nicholls et. al. (PRB 99.094105, 2019).                                   !
    !===========================================================================!

    !===========================================================================!
    subroutine phonon_eels_calculate_thermal_occupations
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                         !
    !                                                                           !
    ! This subroutine calculates the thermal occupation of each phonon mode,    !
    ! according to the Bose-Einstein distribution.                              !
    !                                                                           !
    ! Calling this will allocate & populate the 'phonon_occupations' array.     !
    !===========================================================================!
      implicit none

      ! Dummy variables.
      real(kind=dp)    :: emin_cutoff, x
      integer          :: qi, m

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Allocate phonon_occupations
      allocate(phonon_occupations(n_branches,n_qpts))

      ! Even at small T, occupations blow up for very low frequency modes.
      ! So specify a cutoff (in meV); below this, occupations are capped.
      emin_cutoff = 1.0E-6_dp

      do m=1,n_branches
        do qi=1,n_qpts
          ! Check energy of this mode.
          if ( ABS(phonon_eigenvalues(m,qi)) < emin_cutoff ) then
            x = emin_cutoff / (meV_per_K * ADP_temperature)
          else
            x = phonon_eigenvalues(m,qi) / (meV_per_K * ADP_temperature)
          ! Now, prevent possible overflow in exponential (common when T is very small).
          if ( ABS(x) > 600.0_dp ) then
            x = 600.0_dp
          end if
          ! Calculate occupation of this mode, using the Bose-Einstein function.
          phonon_occupations(m,qi) = 1 / ( EXP(x) - 1 )
          end if
        end do
      end do

      calculated_thermal_occupations = .true.

    end subroutine phonon_eels_calculate_thermal_occupations

    !===========================================================================!
    subroutine phonon_eels_impact_calculate_intensity
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine calculates the EELS intensity spectrum, as given by       !
    ! Eqn. 10 in Nicholls et. al. (PRB 99.094105, 2019).                        !
    !                                                                           !
    ! Note that in the paper, there should be a delta-function in Eqn. 10       !
    ! that enforces q0 + q = RLV. Since we're restricting to the first BZ, this !
    ! effectively translates into q0 = q everywhere (and that is what is        !
    ! calculated here).                                                         !
    !                                                                           !
    ! The inverse square root mass factor is not included here, because we use  !
    ! mass-unweighted eigenvectors in this module (equivalently; the inverse    !
    ! square root mass factor has already been applied at the parsing stage).   !
    !                                                                           !
    ! Calling this will allocate & populate the 'eels_intensity' array.         !
    !===========================================================================!
      implicit none

      ! Dummy variables.
      real(kind=dp)    :: qa, qb, qc, qx, qy, qz, q_nrm, aff, dwf
      integer          :: qi, m, ai, Z, eels_intensity_out_file_unit
      complex(kind=dp) :: f, g, h, intensity

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Allocate eels_intensity.
      allocate(eels_intensity(n_branches,n_qpts))

      ! Get thermal occupations.
      CALL phonon_eels_calculate_thermal_occupations

      ! Loop over q-points and branches.
      do qi=1,n_qpts
        ! Get the Cartesian representation of this q-point; will use for dot-product momentarily.
        qa = qpoint_positions(1,qi)
        qb = qpoint_positions(2,qi)
        qc = qpoint_positions(3,qi)
        qx = qa * recip_lattice(1,1) + qb * recip_lattice(2,1) + qc * recip_lattice(3,1)
        qy = qa * recip_lattice(1,2) + qb * recip_lattice(2,2) + qc * recip_lattice(3,2)
        qz = qa * recip_lattice(1,3) + qb * recip_lattice(2,3) + qc * recip_lattice(3,3)
        q_nrm = SQRT(qx**2 + qy**2 + qz**2)
        do m=1,n_branches
          intensity = cmplx_0
          ! Avoid divergence at q=0.
          if (q_nrm < 1.0E-6_dp) then
            eels_intensity(m,qi) = intensity
            cycle
          end if
          do ai=1,n_ions
            ! First bracket.
            CALL get_atomic_number(atomic_species(ai),Z)
            CALL phonon_eels_calculate_aff(ai,qpoint_positions(:,qi),aff)
            CALL phonon_eels_calculate_dwf(ai,qpoint_positions(:,qi),dwf)
            f = ( REAL(Z,kind=dp) - aff * ( REAL(Z,kind=dp) - partial_charges(ai) ) / REAL(Z,kind=dp) ) * dwf
            ! q.eigvec (this is done in Cartesian coordinates).
            g = qx * phonon_eigenvectors(m,ai,1,qi) + qy * phonon_eigenvectors(m,ai,2,qi) + qz * phonon_eigenvectors(m,ai,3,qi)
            ! q.posn (this is done in fractional coordinates).
            h = twopi * ( qa * atomic_positions(1,ai) + qb * atomic_positions(2,ai) + qc * atomic_positions(3,ai) )
            ! Bring it all together.
            intensity = intensity + f * g * EXP(cmplx_i * h)
          end do
          ! Phonon freqs should always be positive, but divide by ABS(freq) just in case.
          !   --> ALTERNATIVE: Set the intensity of negative-frequency modes to zero.
          intensity = intensity * CONJG(intensity) / ( ABS(phonon_eigenvalues(m,qi)) * q_nrm**4 )
          eels_intensity(m,qi) = intensity * (phonon_occupations(m,qi) + 1)
        end do
      end do

      ! Normalize.
      eels_intensity = eels_intensity / MAXVAL(ABS(eels_intensity))

      ! Write to file.
      if (write_eels_intensity_to_file) then
        eels_intensity_out_file_unit = io_file_unit()
        open(newunit=eels_intensity_out_file_unit,file=TRIM(seedname)//"_impact-intensity.dat",status="replace", &
        action="write",form="formatted")
        write(eels_intensity_out_file_unit,*) ""
        write(eels_intensity_out_file_unit,*) " Impact-EELS Intensity"
        write(eels_intensity_out_file_unit,*) ""
        write(eels_intensity_out_file_unit,*) " All energies in meV"
        write(eels_intensity_out_file_unit,*) ""
        CALL phonon_eels_write_hsps_header(eels_intensity_out_file_unit)
        write(eels_intensity_out_file_unit,*) ""
        write(eels_intensity_out_file_unit,*) ""
        write(eels_intensity_out_file_unit,'(A)') "  q_idx    mode_idx      Energy               I_eels"
        write(eels_intensity_out_file_unit,*) ""
        do qi=1,n_qpts
          do m=1,n_branches
            write(eels_intensity_out_file_unit,'(I5,A,I3,A,F12.6,A,E18.12)') qi,"       ",m,"     ",phonon_eigenvalues(m,qi), &
            "        ",REAL(eels_intensity(m,qi), kind=dp)
          end do
        end do
        close(eels_intensity_out_file_unit)
      end if

      calculated_eels_intensity = .true.

    end subroutine phonon_eels_impact_calculate_intensity

    !==================================================================================!
    subroutine phonon_eels_impact_calculate_heatmap
    !==================================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                                !
    !                                                                                  !
    ! This calculates a heatmap-style phonon dispersion for impact vib-EELS, by        !
    ! placing a (Gaussian approximated) delta-function at each (q,j).                  !
    !                                                                                  !
    ! Calling this will allocate & populate the 'impact_heatmap' array.                !
    !==================================================================================!
      implicit none

      ! Dummy variables.
      integer                      :: p, qi, m, heatmap_out_file_unit
      real(kind=dp)                :: delta, f, nrm, x

      ! Allocate impact_heatmap
      allocate(impact_heatmap(n_qpts,N_energies))

      ! This is a prerequisite for calculating the impact heatmap.
      if (.not. calculated_eels_intensity) then
        CALL phonon_eels_impact_calculate_intensity
      end if

      ! Compute impact heatmap.
      do qi=1,n_qpts
        do p=1,N_energies
          f = 0.0_dp
          do m=1,n_branches
            ! Argument of delta-functional exponential.
            x = ( (energies(p) - phonon_eigenvalues(m,qi)) / vibeels_impact_broadening )**2
            ! Prevent possible underflow from exponential result.
            if (x > 500.0_dp) then
                delta = 0.0_dp
            else
                delta = 1 / (vibeels_impact_broadening * SQRT(pi)) * EXP(-x)
            end if
            f = f + eels_intensity(m,qi) * delta
            impact_heatmap(qi,p) = f
          end do
        end do
      end do

      ! Normalize.
      nrm = MAXVAL(impact_heatmap)
      impact_heatmap = impact_heatmap / nrm

      ! Write heatmap to file.
      if (write_impact_heatmap_to_file) then
        heatmap_out_file_unit = io_file_unit()
        open(newunit=heatmap_out_file_unit,file=TRIM(seedname)//"_impact-heatmap.dat",status="replace", &
        action="write",form="formatted")
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,'(A)') " Heatmap for Impact vib-EELS"
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,'(A,F6.3)') " Broadening (meV): ",REAL(vibeels_impact_broadening, kind=dp)
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,'(A)') " All energies in meV"
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,'(A,I5)') " Number of q-points                 : ",n_qpts
        write(heatmap_out_file_unit,'(A,I5)') " Number of energies at each q-point : ",N_energies
        write(heatmap_out_file_unit,'(A)') ""
        CALL phonon_eels_write_hsps_header(heatmap_out_file_unit)
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,'(A)') ""
        write(heatmap_out_file_unit,"(A,A,A,A,A)") "   q_idx","             ","Energy","                 ","Total Loss"
        write(heatmap_out_file_unit,'(A)') ""
        do qi=1,n_qpts
          do p=1,N_energies
            write(heatmap_out_file_unit,'(A,I5,A,ES16.8,A,ES16.8E3)') " ",qi,"         ",energies(p), &
            "         ",impact_heatmap(qi,p)
          end do
          write(heatmap_out_file_unit,"(A)")
        end do
        close(heatmap_out_file_unit)
      end if

      calculated_impact_heatmap = .true.

    end subroutine phonon_eels_impact_calculate_heatmap

    !===========================================================================!
    subroutine phonon_eels_calculate_dwf(atom_idx,qpt_frac,dwf)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine calculates the Debye-Waller factor, specifically the      !
    ! amplitude EXP(-W), for 'atom_idx' at the q-vector 'qpt_frac'. The result  !
    ! is stored in 'dwf'.                                                       !
    !===========================================================================!
      implicit none

      ! Arguments.
      integer, INTENT(IN)        :: atom_idx
      real(kind=dp), INTENT(IN)  :: qpt_frac(3)
      real(kind=dp), INTENT(OUT) :: dwf

      ! Dummy variables.
      integer          :: i, j
      real(kind=dp)    :: W
      real(kind=dp)    :: q_cart(3)

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! q-vector in Cartesian basis.
      q_cart(1) = qpt_frac(1) * recip_lattice(1,1) + qpt_frac(2) * recip_lattice(2,1) + qpt_frac(3) * recip_lattice(3,1)
      q_cart(2) = qpt_frac(1) * recip_lattice(1,2) + qpt_frac(2) * recip_lattice(2,2) + qpt_frac(3) * recip_lattice(3,2)
      q_cart(3) = qpt_frac(1) * recip_lattice(1,3) + qpt_frac(2) * recip_lattice(2,3) + qpt_frac(3) * recip_lattice(3,3)

      ! Compute DWF.
      W = 0.0_dp
      do i=1,3
        do j=1,3
          W = W + 0.5_dp * q_cart(i) *  atomic_displacement_params(atom_idx,i,j) * q_cart(j)
        end do
      end do
      dwf = EXP(-W)

    end subroutine phonon_eels_calculate_dwf

    !===========================================================================!
    subroutine phonon_eels_calculate_aff(atom_idx,qpt_frac,aff)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine calculates the atomic form factor (AFF) for 'atom_idx'    !
    ! at the q-vector 'qpt_frac'. The result is stored in 'aff'.                !
    !                                                                           !
    ! By default, X-ray AFFs are calculated from the Cromer-Mann coefficients.  !
    ! If electronic AFFs have instead been requested by the user, then          !
    ! the electronic scattering factor is first calculated using the Peng       !
    ! coefficients (L. M. Peng, 1999), and then converted back to an AFF via    !
    ! the Mott-Bethe relation.                                                  !
    !                                                                           !
    ! The specific analytic form of the Mott-Bethe relation as implemented      !
    ! here is taken from:                                                       !
    !                                                                           !
    !                 'Reflection High-Energy Electron Diffraction'             !
    !                     A. Ichimiya and P. I. Cohen (2004)                    !
    !                         Ch. 9, Page 117, Eqn. 9.32                        !
    !===========================================================================!
      implicit none

      ! Subroutine arguments.
      integer, INTENT(IN)        :: atom_idx
      real(kind=dp), INTENT(IN)  :: qpt_frac(3)
      real(kind=dp), INTENT(OUT) :: aff

      ! Dummy variables.
      real(kind=dp) :: a(4), b(4), c, qx, qy, qz, q_nrm, s
      real(kind=dp) :: x(4), electronic_aff_prefactor
      integer       :: Z, i

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Compute aff.
      CALL get_atomic_number(atomic_species(atom_idx),Z)
      if (vibeels_use_electronic_affs) then
        CALL get_peng_electron_coeff(Z,a,b,c)
      else
        call get_cromer_mann_coeff(Z,a,b,c)
      end if

      ! Get the magnitude of this q-vector, in units of inverse Angstrom.
      ! This is required for use with the fitted parametrization from the tables.
      qx = qpt_frac(1) * recip_lattice(1,1) + qpt_frac(2) * recip_lattice(2,1) + qpt_frac(3) * recip_lattice(3,1)
      qy = qpt_frac(1) * recip_lattice(1,2) + qpt_frac(2) * recip_lattice(2,2) + qpt_frac(3) * recip_lattice(3,2)
      qz = qpt_frac(1) * recip_lattice(1,3) + qpt_frac(2) * recip_lattice(2,3) + qpt_frac(3) * recip_lattice(3,3)
      q_nrm = SQRT( qx**2 + qy**2 + qz**2 )
      s = q_nrm / (4.0_dp*pi)
      aff = 0.0_dp
      do i = 1,4
        aff = aff + a(i) * exp( -b(i) * s**2 )
      end do
      aff = aff + c

      ! If electronic AFFs have been requested, then we must apply the Mott-Bethe relation to obtain an AFF.
      if (vibeels_use_electronic_affs) then
        electronic_aff_prefactor = 0.023934_dp
        aff = REAL(Z,kind=dp) - aff * s**2 / electronic_aff_prefactor
      end if

    end subroutine phonon_eels_calculate_aff

    !==================================================================================!
    subroutine phonon_eels_write_affs_and_dwfs_to_file
    !==================================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                              !
    !                                                                                  !
    ! This subroutine calculates the Atomic Form Factors (AFFs) and Debye-Waller       !
    ! Factors (DWFs) for all atoms, at all of the q-points in the .phonon file, and    !
    ! writes them to file.                                                             !
    !                                                                                  !
    ! Calling this will create 'seedname_aff.dat' and 'seedname_dwf.dat' files.        !
    !==================================================================================!
      implicit none

      ! Dummy variables.
      real(kind=dp)                        :: aff,dwf
      integer                              :: ai,qi,aff_out_file_unit,dwf_out_file_unit
      real(kind=dp), allocatable, save     :: affs_all(:,:)                               ! (atom_idx,qi)
      real(kind=dp), allocatable, save     :: dwfs_all(:,:)                               ! (atom_idx,qi)
      character(len=10)                    :: istr

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Allocate affs_all and dwfs_all.
      allocate(affs_all(n_ions,n_qpts))
      allocate(dwfs_all(n_ions,n_qpts))

      ! Loop over all atoms and q-vectors in the .phonon file.
      do ai=1,n_ions
        do qi=1,n_qpts
          CALL phonon_eels_calculate_aff(ai,qpoint_positions(:,qi),aff)
          affs_all(ai,qi) = aff
          CALL phonon_eels_calculate_dwf(ai,qpoint_positions(:,qi),dwf)
          dwfs_all(ai,qi) = dwf
        end do
      end do

      ! Write AFFs to file.
      ! There one set of AFF data for each element type.
      aff_out_file_unit = io_file_unit()
      open(newunit=aff_out_file_unit,file=TRIM(seedname)//"_aff.dat",status="replace",action="write",form="formatted")
      write(aff_out_file_unit,'(A)') ""
      write(aff_out_file_unit,'(A)') " Atomic Form Factors (AFFs)"
      write(aff_out_file_unit,'(A)') ""
      if (vibeels_use_electronic_affs) then
        write(aff_out_file_unit,'(A)') " AFF Type : Electronic"
      else
        write(aff_out_file_unit,'(A)') " AFF Type : X-Ray"
      end if
      write(aff_out_file_unit,'(A)') ""
      CALL phonon_eels_write_hsps_header(aff_out_file_unit)
      write(aff_out_file_unit,'(A)') ""
      write(aff_out_file_unit,'(A)') ""
      write(aff_out_file_unit,*) "type     q_idx       qa           qb           qc            aff(|q|)"
      write(aff_out_file_unit,'(A)') ""
      do ai=1,n_ions
        ! Skip duplicate elements.
        if (any(atomic_species(:ai-1) == atomic_species(ai))) then
          cycle
        end if
        do qi=1,n_qpts
          write(aff_out_file_unit,'(A,A,I5,A,F8.5,A,F8.5,A,F8.5,A,ES16.8)') " "//atomic_species(ai),"     ",qi,"     ", &
          qpoint_positions(1,qi),"     ",qpoint_positions(2,qi),"     ",qpoint_positions(3,qi),"     ",affs_all(ai,qi)
        end do
        write(aff_out_file_unit,*) ""
      end do
      close (aff_out_file_unit)

      ! Write DWFs to file.
      ! There is one set of DWF data for each atom in the unit cell.
      dwf_out_file_unit = io_file_unit()
        open(newunit=dwf_out_file_unit,file=TRIM(seedname)//"_dwf.dat",status="replace",action="write",form="formatted")
        write(dwf_out_file_unit,'(A)') ""
        write(dwf_out_file_unit,'(A)') " Debye-Waller Factors (DWFs)"
        write(aff_out_file_unit,'(A)') ""
        CALL phonon_eels_write_hsps_header(dwf_out_file_unit)
        write(dwf_out_file_unit,'(A)') ""
        write(aff_out_file_unit,'(A)') ""
        write(dwf_out_file_unit,*) "atom     q_idx      qa           qb           qc              dwf(q)"
        write(dwf_out_file_unit,'(A)') ""
        do ai=1,n_ions
          write(istr,'(I3)') ai
          do qi=1,n_qpts
            write(dwf_out_file_unit,'(A,A,I5,A,F8.5,A,F8.5,A,F8.5,A,ES16.8)') " "//trim(adjustl(atomic_species(ai)))// &
            trim(adjustl(istr)),"     ",qi,"     ",qpoint_positions(1,qi),"     ",qpoint_positions(2,qi), "     ", &
            qpoint_positions(3,qi),"     ",dwfs_all(ai,qi)
          end do
          write(dwf_out_file_unit,*) ""
        end do
      close (dwf_out_file_unit)

      ! Deallocate temporary arrays.
      deallocate(affs_all)
      deallocate(dwfs_all)

    end subroutine phonon_eels_write_affs_and_dwfs_to_file








    !==================================================================================!
    !                              *** ALOOF ROUTINES ***                              !
    !                                                                                  !
    ! The routines here are used for the aloof method, as described in Radtke et. al.  !
    ! (PRL 119.027402, 2017).                                                          !
    !==================================================================================!

    !==================================================================================!
    subroutine phonon_eels_aloof_calculate_loss_probability
    !==================================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                              !
    !                                                                                  !
    ! This subroutine calculates the aloof loss probability, according to Eqn. 3 of    !
    ! Radtke et. al. (PRL 119.027402, 2017).                                           !
    !                                                                                  !
    ! At the moment, the loss probability is normalized such that its maximum is 1,    !
    ! and this is what is written to the output file.                                  !
    !                                                                                  !
    ! Calling this will allocate & populate the 'aloof_loss_probability' array.        !
    !==================================================================================!
      implicit none

      ! Dummy variables.
      integer                      :: w, wc, p, aloof_loss_probability_out_file_unit
      real(kind=dp)                :: beam_energy_au, rel_beta, rel_gamma, electron_velocity_au, impact_parameter_au
      real(kind=dp)                :: kx, ky, K, dphi, prefactor, f, lor_broadening_au, f_lor, dE

      ! Temporary array, used to broaden the loss.
      real(kind=dp),allocatable :: aloof_loss_broadened(:)

      ! This is a prerequisite for calculating the aloof loss probability.
      if (.not. calculated_polarizability) then
        CALL phonon_eels_aloof_calculate_polarizability
      end if

      ! Allocate aloof_loss_probability.
      allocate(aloof_loss_probability(N_energies))

      ! Allocate aloof_loss_broadened (temp array), which is used to broaden the loss.
      allocate(aloof_loss_broadened(N_energies))

      ! Get dphi (integration weight for angular integration).
      dphi = vibeels_aloof_phi_spacing

      ! Convert electron beam energy in units of keV, to electron velocity in Hartree atomic units.
      beam_energy_au = vibeels_aloof_electron_beam_energy * 1.0E3_dp * eV_to_hartree
      rel_gamma = 1.0_dp + beam_energy_au * fine_structure_constant**2
      rel_beta = sqrt(1.0_dp - 1.0_dp / rel_gamma**2)
      electron_velocity_au = rel_beta / fine_structure_constant

      ! Convert impact parameter from nanometre to Bohr radii.
      impact_parameter_au = vibeels_aloof_impact_parameter * nanometre2bohr

      ! Prefactor for aloof loss expression.
      ! We normalize to 1 later, so this doesn't do anything useful at the moment.
      prefactor = 1.00 / electron_velocity_au**2

      do w=1,N_energies
        f = 0.0_dp
        do p=1,N_phi_values
          kx = omegas_hartree(w) / electron_velocity_au
          ky = kx * TAN(phis(p))
          K = SQRT(kx**2 + ky**2)
          f = f + AIMAG(polarizability(w,p)) * EXP(-2.0_dp * impact_parameter_au * K) / COS(phis(p)) * dphi
        end do
        aloof_loss_probability(w) = f * prefactor
      end do

      ! Now, we convolve the aloof loss spectrum with a Lorentzian.

      ! Convert broadening to au.
      lor_broadening_au = vibeels_aloof_loss_broadening * meV_to_hartree

      ! Perform the convolution.
      dE = omegas_hartree(2) - omegas_hartree(1)
      do w=1,N_energies
        aloof_loss_broadened(w) = 0.0_dp
        do wc=1,N_energies
          f_lor = (lor_broadening_au/2) / ( pi * ( omegas_hartree(w)-omegas_hartree(wc) )**2 + (lor_broadening_au/2)**2 )
          aloof_loss_broadened(w) = aloof_loss_broadened(w) + aloof_loss_probability(wc) * f_lor * dE
        end do
      end do

      ! Update the loss spectrum with its broadened version.
      aloof_loss_probability = aloof_loss_broadened

      ! We can deallocate aloof_loss_broadened now.
      deallocate(aloof_loss_broadened)

      ! Normalize the loss to 1.
      aloof_loss_probability = aloof_loss_probability / MAXVAL(ABS(aloof_loss_probability))

      ! Write to file.
      if (write_aloof_loss_probability_to_file) then
        aloof_loss_probability_out_file_unit = io_file_unit()
        open(newunit=aloof_loss_probability_out_file_unit,file=TRIM(seedname)//"_aloof-loss.dat",status="replace", &
        action="write",form="formatted")
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A)') " Aloof-EELS Loss Probability"
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,A,A)') " Crystal Surface : ",vibeels_aloof_surface_plane,"-plane"
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,A,A)') " Beam Direction  : ",vibeels_aloof_beam_direction,"-axis"
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,F12.6)') " Electron Beam Energy (keV)            : ", &
        vibeels_aloof_electron_beam_energy
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,F12.6)') " Aloof Impact Parameter (nanometre)    : ", &
        vibeels_aloof_impact_parameter
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,F12.6)') " Intrinsic Broadening (meV)            : ", &
        vibeels_aloof_intrinsic_broadening
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A,F12.6)') " Aloof Loss Spectrum Broadening (meV)  : ", &
        vibeels_aloof_loss_broadening
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A)') " All energies in meV"
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        write(aloof_loss_probability_out_file_unit,'(A)') "    Energy                 P_loss"
        write(aloof_loss_probability_out_file_unit,'(A)') ""
        do w=1,N_energies
          write(aloof_loss_probability_out_file_unit,'(F12.6,A,E18.12)') energies(w),"         ", &
          REAL(aloof_loss_probability(w), kind=dp)
        end do
        close(aloof_loss_probability_out_file_unit)
      end if

      calculated_aloof_loss_probability = .true.

    end subroutine phonon_eels_aloof_calculate_loss_probability

    !==================================================================================!
    subroutine phonon_eels_aloof_calculate_polarizability
    !==================================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                              !
    !                                                                                  !
    ! This subroutine calculates the polarizability, as given by Eqn. 2 in Radtke      !
    ! et. al. (PRL 119.027402, 2017).                                                  !
    !                                                                                  !
    ! The original work of Radtke et. al. assumes a diagonal eps_lf; if                !
    ! 'vibeels_aloof_include_offdiag' is .true., then the fully-general case is        !
    ! instead treated, and the off-diagonal components of eps_lf are included.         !
    !                                                                                  !
    ! Calling this will allocate & populate the 'polarizability' array.                !
    !==================================================================================!
      implicit none

      ! Dummy variables.
      integer           :: i, j, w, p, alpha_out_file_unit
      real(kind=dp)     :: kx,ky
      complex(kind=dp)  :: f,n_eff,n_eff_prev

      ! This is a prerequisite for calculating the polarizability.
      if (.not. calculated_lf_eps) then
        CALL phonon_eels_calculate_lf_dielectric_tensor
      end if

      ! Allocate polarizability array.
      allocate(polarizability(N_energies,N_phi_values))

      do w=1,N_energies
        do p=1,N_phi_values
          ! Relative components of wavevector.
          kx = COS(phis(p))
          ky = SIN(phis(p))
          ! Now, compute alpha.
          !**********************************************************************************************************************
          ! Diagonal case (original, from Radtke).
          if (.not. vibeels_aloof_include_offdiag) then
            f = lf_dielectric_tensor(3,3,w) * ( kx**2 * lf_dielectric_tensor(1,1,w) + ky**2 * lf_dielectric_tensor(2,2,w) )
          ! Fully general case, that allows for non-diagonal eps_lf.
          else
            f = lf_dielectric_tensor(3,3,w) * ( kx**2 * lf_dielectric_tensor(1,1,w) + 2*kx*ky * lf_dielectric_tensor(1,2,w) &
            + ky**2 * lf_dielectric_tensor(2,2,w) ) - ( kx * lf_dielectric_tensor(1,3,w) + ky * lf_dielectric_tensor(2,3,w) )**2
          end if
          !**********************************************************************************************************************
          n_eff = SQRT(f)
          ! Optional branch cut check; fixes discontinuities in alpha.
          if (alpha_enforce_analytic_continuity) then
            if (w > 1) then
              if (abs(n_eff - n_eff_prev) > abs(-n_eff - n_eff_prev)) then
                  n_eff = -n_eff
              end if
            end if
            n_eff_prev = n_eff
          end if
          ! Finally, calculate alpha.
          polarizability(w,p) = (n_eff - 1) / (n_eff + 1)
        end do
      end do

      ! Write polarizability to file.
      if (write_alpha_lf_to_file) then
        alpha_out_file_unit = io_file_unit()
        open(newunit=alpha_out_file_unit,file=TRIM(seedname)//"_alpha.dat",status="replace",action="write",form="formatted")
        write(alpha_out_file_unit,*) ""
        write(alpha_out_file_unit,*) "Polarizability (alpha)"
        write(alpha_out_file_unit,'(A)') ""
        write(alpha_out_file_unit,'(A,A,A)') " Crystal Surface : ",vibeels_aloof_surface_plane,"-plane"
        write(alpha_out_file_unit,'(A)') ""
        write(alpha_out_file_unit,'(A,A,A)') " Beam Direction  : ",vibeels_aloof_beam_direction,"-axis"
        write(alpha_out_file_unit,'(A)') ""
        write(alpha_out_file_unit,"(A)") " Phi = Angle between kx and ky (in radians)"
        write(alpha_out_file_unit,*) ""
        write(alpha_out_file_unit,'(A)') " All energies in meV"
        write(alpha_out_file_unit,*) ""
        write(alpha_out_file_unit,'(A)') ""
        write(alpha_out_file_unit,'(A)') "    Phi                 Energy                 Real(alpha)                Imag(alpha)"
        write(alpha_out_file_unit,*) ""
        do p=1,N_phi_values
          do w=1,N_energies
            write(alpha_out_file_unit,'(F12.6,A,F12.6,A,E18.12,A,E18.12)') phis(p),"        ", &
            omegas_hartree(w)/meV_to_hartree, "           ", REAL(polarizability(w,p), kind=dp),"         ", &
            AIMAG(polarizability(w,p))
          end do
          write(alpha_out_file_unit,*) ""
        end do
        close(alpha_out_file_unit)
      end if

      calculated_polarizability = .true.

    end subroutine phonon_eels_aloof_calculate_polarizability

    !===========================================================================!
    subroutine phonon_eels_calc_oscillator_strength_tensor
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine calculates the oscillator strength tensor, as given by    !
    ! Eqn. 54 of Gonze & Lee (PRB 55.10355, 1997).                              !
    !                                                                           !
    ! As per Gonze & Lee, mass-unweighted eigenvectors are assumed.             !
    !                                                                           !
    ! Calling this will allocate & populate the 'osc_strength_tensor' array.    !
    !===========================================================================!
      implicit none
    
      ! Dummy variables.
      integer                       :: m, a, b, k, ap, kp, bp, osc_strength_tensor_out_file_unit
      complex(kind=dp)              :: f

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Allocate osc_strength_tensor.
      allocate(osc_strength_tensor(n_branches,3,3))

      ! Calculate oscillator strength tensor.
      do m=1,n_branches
        do a=1,3
          do b=1,3
            f = cmplx_0
            do k=1,n_ions
              do ap=1,3
                do kp=1,n_ions
                  do bp=1,3
                    f = f + born_eff_ch_tensor(k,a,ap) * CONJG(gamma_eigenvectors(m,k,ap)) &
                    * born_eff_ch_tensor(kp,b,bp) * gamma_eigenvectors(m,kp,bp)
                  end do
                end do
              end do
            end do
            osc_strength_tensor(m,a,b) = f
          end do
        end do
      end do

      ! Write to file.
      if (write_osc_strength_tensor_to_file) then
        osc_strength_tensor_out_file_unit = io_file_unit()
        open(newunit=osc_strength_tensor_out_file_unit,file=TRIM(seedname)//"_ost.dat",status="replace", &
        action="write",form="formatted")
        write(osc_strength_tensor_out_file_unit,'(A)') ""
        write(osc_strength_tensor_out_file_unit,'(A)') " Oscillator Strength Tensor (at Gamma-point)"
        write(osc_strength_tensor_out_file_unit,'(A)') ""
        write(osc_strength_tensor_out_file_unit,'(A)') " --> Calculated using Eqn. 54 of Gonze & Lee (PRB 55.10355, 1997)"
        write(osc_strength_tensor_out_file_unit,'(A)') ""
        do m=1,n_branches
          write(osc_strength_tensor_out_file_unit,"(A,I3,A,F16.10,A)") "  Oscillator Strength Tensor for Mode ",m, &
          ", frequency: ",gamma_eigenvalues(m), " meV"
          write(osc_strength_tensor_out_file_unit,"(ES14.6,ES14.6,ES14.6)") REAL(osc_strength_tensor(m,1,1)), &
          REAL(osc_strength_tensor(m,1,2)), REAL(osc_strength_tensor(m,1,3))
          write(osc_strength_tensor_out_file_unit,"(ES14.6,ES14.6,ES14.6)") REAL(osc_strength_tensor(m,2,1)), &
          REAL(osc_strength_tensor(m,2,2)), REAL(osc_strength_tensor(m,2,3))
          write(osc_strength_tensor_out_file_unit,"(ES14.6,ES14.6,ES14.6)") REAL(osc_strength_tensor(m,3,1)), &
          REAL(osc_strength_tensor(m,3,2)), REAL(osc_strength_tensor(m,3,3))
          write(osc_strength_tensor_out_file_unit,'(A)') ""
        end do
        close(osc_strength_tensor_out_file_unit)
      end if

      calculated_osc_strength_tensor = .true.

    end subroutine phonon_eels_calc_oscillator_strength_tensor

    !===========================================================================!
    subroutine phonon_eels_calc_mode_resolved_born_tensor
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine calculates the mode-resolved Born effective charge        !
    ! tensor, as given by Eqn. 53 of Gonze & Lee (PRB 55.10355, 1997).          !
    !                                                                           !
    ! As per Gonze & Lee, mass-unweighted eigenvectors are assumed.             !
    !                                                                           !
    ! This subroutine is not used at present for the calculation of any aloof   !
    ! vib-EELS quantities (rather, the oscillator strength tensor is used       !
    ! instead).                                                                 !
    !                                                                           !
    ! Calling this will allocate & populate the                                 !
    ! 'mode_resolved_born_eff_ch_tensor' array.                                 !
    !===========================================================================!
      implicit none

      ! Dummy variables.
      integer           :: m, k, a, b
      complex(kind=dp)  :: f, norm

      ! Check that we have all required quantities.
      if (.not. phonon_eels_prep_done) then
        CALL phonon_eels_prepare
      end if

      ! Allocate mode_resolved_born_eff_ch_tensor.
      allocate(mode_resolved_born_eff_ch_tensor(n_branches,3))

      ! Now, populate mode_resolved_born_eff_ch_tensor.
      do m=1,n_branches
        do a=1,3
          f = cmplx_0
          norm = cmplx_0
          do k=1,n_ions
            do b=1,3
              f = f + born_eff_ch_tensor(k,a,b) * gamma_eigenvectors(m,k,b)
              norm = norm + CONJG( gamma_eigenvectors(m,k,b) ) * gamma_eigenvectors(m,k,b)
            end do
          end do
          mode_resolved_born_eff_ch_tensor(m,a) = f / SQRT(norm)
        end do
      end do

      ! For debug: Print mode-resolved BEC tensor (based on iprint level).
      if (iprint > 1) then
        write(stdout,*) ""
        write(stdout,*) "MODE-RESOLVED BEC TENSOR (at Gamma-point):"
        write(stdout,*) ""
        do m = 1, n_branches
          write(stdout,*) "Resolved BEC Tensor for Mode ", m
          write(stdout,"(3('(',ES14.6,',',ES14.6,')'))") &
            ( REAL(mode_resolved_born_eff_ch_tensor(m,a), kind=dp), &
              AIMAG(mode_resolved_born_eff_ch_tensor(m,a)), a=1,3 )
          write(stdout,*) ""
        end do
      end if

      calculated_mode_resolved_bec_tensor = .true.

    end subroutine phonon_eels_calc_mode_resolved_born_tensor

    !==================================================================================!
    subroutine phonon_eels_calculate_lf_dielectric_tensor
    !==================================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                              !
    !                                                                                  !
    ! This subroutine calculates the low-frequency, omega-dependent dielectric         !
    ! tensor, as given by Eqn. 55 of Gonze & Lee (PRB 55.10355, 1997).                 !
    !                                                                                  !
    ! Calling this will allocate & populate the 'lf_dielectric_tensor' array.          !
    !                                                                                  !
    ! NOTE 1: Hartree atomic units are used internally in this subroutine, because     !
    !         Eqn. 55 in Gonze & Lee requires it.                                      !
    !                                                                                  !
    ! NOTE2 : For now, the intrinsic broadening (supplied in .odi file) is a single    !
    !         real value, i.e. the same for all modes. In general, this should be      !
    !         mode-specific, and can be calculated with an anharmonic calculation.     !
    !==================================================================================!
      implicit none

      ! Dummy variables.
      integer                      :: mode, atom, i, j, p
      integer                      :: eps_lf_out_file_unit, first_cpt, second_cpt
      character(len=1)             :: cpts_string(3)
      real(kind=dp)                :: broadening_au, cell_volume_au
      real(kind=dp),allocatable    :: gamma_eigvals_hartree(:)
      complex(kind=dp)             :: eps_lf, numerator, denominator

      ! This is a prerequisite for calculating eps_lf.
      if (.not. calculated_osc_strength_tensor) then
        CALL phonon_eels_calc_oscillator_strength_tensor
      end if

      ! Convert intrinsic broadening (linewidth) to Hartree.
      broadening_au = vibeels_aloof_intrinsic_broadening * meV_to_hartree

      ! Convert cell volume to Bohr^3.
      cell_volume_au = cell_volume * ang2bohr**3

      ! Gamma-point phonon frequencies in Hartree.
      allocate(gamma_eigvals_hartree(n_branches))
      gamma_eigvals_hartree = gamma_eigenvalues * meV_to_hartree

      ! Allocate lf_dielectric_tensor.
      allocate(lf_dielectric_tensor(3,3,N_energies))

      ! Calculate low-frequency dielectric tensor.
      do i=1,3
        do j=1,3
          do p=1,N_energies
            eps_lf = inf_dielectric_tensor(i,j)
            do mode=1,n_branches
              if ( eps_lf_skip_acoustic_modes .and. (mode .lt. 4) ) cycle ! Optionally skip acoustic modes
              numerator = osc_strength_tensor(mode,i,j)
              denominator = gamma_eigvals_hartree(mode)**2 - (omegas_hartree(p) + cmplx_i*broadening_au)**2
              eps_lf = eps_lf + (fourpi/cell_volume_au) * numerator/denominator
            end do
            lf_dielectric_tensor(i,j,p) = eps_lf
          end do
        end do
      end do

      ! Deallocate temporary arrays.
      deallocate(gamma_eigvals_hartree)

      ! Write eps_lf to file.
      ! This will write out a separate file for each of the 'xx', 'xy', 'xz', 'yy', 'yz', 'zz' components.
      if (write_eps_lf_to_file) then
        cpts_string(1) = "x"
        cpts_string(2) = "y"
        cpts_string(3) = "z"
        eps_lf_out_file_unit = io_file_unit()
        open(newunit=eps_lf_out_file_unit,file=TRIM(seedname)//"_eps-lf.dat", status="replace",action="write",form="formatted")
        write(eps_lf_out_file_unit,'(A)') ""
        write(eps_lf_out_file_unit,'(A)') " Low-frequency dielectric function"
        if (aloof_rotation_applied) then
          write(eps_lf_out_file_unit,'(A)') ""
          write(eps_lf_out_file_unit,'(A)') " NOTE: For this aloof calculation, a rotation was applied to the unit cell"
          write(eps_lf_out_file_unit,'(A,A,A)') "       --> The normal to the ",vibeels_aloof_surface_plane, &
          "-plane is now parallel to the z-direction"
          write(eps_lf_out_file_unit,'(A,A,A)') "       --> The crystal ",vibeels_aloof_beam_direction, &
          "-axis is now parallel to the x-direction"
        end if
        write(eps_lf_out_file_unit,'(A)') ""
        write(eps_lf_out_file_unit,'(A, F6.3)') " Intrinsic Broadening (meV) : ",vibeels_aloof_intrinsic_broadening
        write(eps_lf_out_file_unit,'(A)') ""
        write(eps_lf_out_file_unit,'(A)') " All energies in meV"
        write(eps_lf_out_file_unit,'(A)') ""
        write(eps_lf_out_file_unit,'(A)') ""
        do second_cpt=1,3
          do first_cpt=1,second_cpt
            write(eps_lf_out_file_unit,'(A,A,A)') "                             Tensor Component : ", &
            cpts_string(first_cpt),cpts_string(second_cpt)
            write(eps_lf_out_file_unit,'(A)') ""
            write(eps_lf_out_file_unit,'(A)') "    Energy                        Real(eps_lf)                    Im(eps_lf)"
            write(eps_lf_out_file_unit,'(A)') ""
            do p=1,N_energies
              write(eps_lf_out_file_unit, '(F12.6, A, ES24.12, A, ES24.12)') omegas_hartree(p)/meV_to_hartree, &
              "            ", REAL(lf_dielectric_tensor(first_cpt,second_cpt,p), kind=dp), "       ", &
              AIMAG(lf_dielectric_tensor(first_cpt,second_cpt,p))
            end do
            write(eps_lf_out_file_unit,'(A)') ""
            write(eps_lf_out_file_unit,'(A)') ""
          end do
        end do
        close(eps_lf_out_file_unit)
      end if

      calculated_lf_eps = .true.

    end subroutine phonon_eels_calculate_lf_dielectric_tensor








    !===========================================================================!
    !                        *** PREPARE EVERYTHING ***                         !
    !===========================================================================!

    !=========================================================================!
    subroutine phonon_eels_prepare
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                       !
    !                                                                         !
    ! This subroutine prepares everything for vib-EELS. It checks that the    !
    ! relevant data has been parsed from the .phonon and .odd files,          !
    ! allocates and sets the frequency arrays, and also sets up the array of  !
    ! phi values, which are used to calculate the polarizability integral for !
    ! the aloof loss probability.                                             !
    !                                                                         !
    ! If so requested, the phonon bands are also appropriately reordered      !
    ! according to eigenvector character, and high-symmetry points along the  !
    ! dispersion path are detected.                                           !
    !                                                                         !
    ! This subroutine should be called before any of the vib-EELS calculation !
    ! subroutines are used (if not, those subroutines will check that it has  !
    ! been called anyway).                                                    !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      integer :: p, m, ai

      ! For rotating the unit cell, if required for aloof.
      real(kind=dp) :: R(3,3), Rf(3,3), Rtot(3,3), identity(3,3)
      real(kind=dp) :: v(3), vx, vy, nrm, theta

      ! Get relevant data from .odd file.
      if (.not. parsed_eels_data_from_odd) then
        CALL odd_read_vib_eels_data
      end if

      ! Parse .phonon file.
      if (.not. parsed_phonon_file) then
        CALL phonon_eels_read_phonon_file
      end if

      ! Reorder phonon bands, if requested.
      if (vibeels_reorder_phonon_bands) then
        CALL phonon_eels_reorder_phonon_bands
      end if

      ! Detect high-symmetry points along dispersion path.
      ! Don't do this for aloof-only mode, as the user may have supplied Gamma-point only data,
      ! so there is no high-symmetry path in the first place.
      if ( index(phonon_eels_task, 'aloof') .eq. 0 ) then
        CALL phonon_eels_detect_hsps
      end if

      ! If vibeels_max_energy was left unset in the .odi file, then it will be negative (by default).
      ! So, we set it here to the highest phonon energy plus 50 %.
      if (vibeels_max_energy < 0.0_dp) then
        vibeels_max_energy = MAXVAL(ABS(phonon_eigenvalues)) * 1.50
      end if

      ! Allocate & populate the energies and omegas_hartree arrays.
      ! This will also set N_energies.
      N_energies = int( (vibeels_max_energy - vibeels_min_energy) / vibeels_spacing )
      allocate(energies(N_energies))
      allocate(omegas_hartree(N_energies))
      do p=1,N_energies
        energies(p) = vibeels_min_energy + (p-1) * vibeels_spacing
        omegas_hartree(p) = energies(p) * meV_to_hartree
      end do

      ! If we're just doing impact vib-EELS, then we're done here.
      if ( index(phonon_eels_task, 'impact') > 0 ) then
        phonon_eels_prep_done = .true.
        return
      end if

      ! Needed for aloof only; allocate phi array, and populate over the range (-pi/2,pi/2)
      ! This will also set N_phi_values.
      N_phi_values = int(pi / vibeels_aloof_phi_spacing)
      allocate(phis(N_phi_values))
      do p=1,N_phi_values
        phis(p) = (REAL(p,kind=dp)/(N_phi_values+1) - 0.5_dp) * pi
      end do

      ! If we're doing any aloof calculations, then we may need to rotate the aloof vector and tensor quantities
      ! appropriately, so that the crystal surface normal is aligned along the z-direction (for the aloof theory
      ! of Radtke et. al., this must be the normal direction).

      ! It is perfectly possible that the cell is already appropriately oriented, in which
      ! case the code below will not modify the original quantities.

      ! First, get the matrix 'R' that rotates the crystal such that the normal to the surface plane is parallel to the z-direction.
      CALL phonon_eels_get_normal_rotation_matrix(vibeels_aloof_surface_plane,R)

      ! We now have the R matrix that rotates the crystal, such that the appropriate RLV is parallel to the z-axis.
      ! Now, we just need to rotate so that the desired real-space axis is parallel to the beam direction.
      ! In the aloof theory of Radtke et. al., this is defined to be the x-direction.
      ! First, get the direction that the desired real-space axis is currently pointing in (after the application of R).
      if (index(vibeels_aloof_beam_direction,"a") > 0) then
        v = matmul(R,real_lattice(1,:))
      else if (index(vibeels_aloof_beam_direction,"b") > 0) then
        v = matmul(R,real_lattice(2,:))
      else if (index(vibeels_aloof_beam_direction,"c") > 0) then
        v = matmul(R,real_lattice(3,:))
      else
        CALL io_error("ERROR: vibeels_aloof_beam_direction must be one of 'a', 'b' or 'c'")
      end if

      ! Get the angle of this vector in the x-y plane.
      theta = ATAN2( v(2),v(1) )

      ! Now construct a rotation about the z-axis that brings this angle back to zero.
      ! That is; rotate by the negative of this angle, around z.
      ! This will then ensure that this axis lies parallel to the x-direction, as desired.
      Rf(1,:) = [cos(theta),sin(theta),0.0_dp]
      Rf(2,:) = [-sin(theta),cos(theta),0.0_dp]
      Rf(3,:) = [0.0_dp,0.0_dp,1.0_dp]

      ! Finally, get the full rotation matrix.
      Rtot = matmul(Rf,R)

      ! Check if Rtot is identity.
      ! If it isn't, then a rotation will be applied to the cell.
      identity = 0.0_dp
      identity(1,1) = 1.0_dp
      identity(2,2) = 1.0_dp
      identity(3,3) = 1.0_dp
      if ( ALL( ABS(Rtot - identity) < 1.0E-4_dp ) ) then
        aloof_rotation_applied = .false.
      else
        aloof_rotation_applied = .true.
      end if

      ! Now, we rotate the quantities required for the aloof theory.
      ! Note that the code below does not touch, or affect, any quantities currently used in the impact theory.
      ! At present therefore, this is a safe approach.

      ! Rotate Gamma-point phonon eigenvectors.
      do m=1,n_branches
        do ai=1,n_ions
          gamma_eigenvectors(m,ai,:) = matmul(Rtot,gamma_eigenvectors(m,ai,:))
        end do
      end do

      ! Rotate inf_dielectric_tensor.
      inf_dielectric_tensor = matmul( Rtot, matmul(inf_dielectric_tensor, TRANSPOSE(Rtot)) )

      ! Rotate born_eff_ch_tensor.
      do ai=1,n_ions
        born_eff_ch_tensor(ai,:,:) = matmul( Rtot, matmul(born_eff_ch_tensor(ai,:,:), TRANSPOSE(Rtot)) )
      end do

      phonon_eels_prep_done = .true.

    end subroutine phonon_eels_prepare

    !=========================================================================!
    subroutine phonon_eels_read_phonon_file
    !=========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                     !
    !                                                                         !
    ! This subroutine reads in a CASTEP 'seedname.phonon' file, and parses    !
    ! the cell vectors, atomic positions, q-points, phonon frequencies, and   !
    ! phonon eigenvectors. The reciprocal lattice vectors are also computed,  !
    ! alongside the volume of the real-space cell.                            !
    !                                                                         !
    ! Calling this will allocate & populate the 'real_lattice',               !
    ! 'recip_lattice', 'atomic_positions', 'qpoint_positions',                !
    ! 'phonon_eigenvalues', 'phonon_eigenvectors', 'gamma_eigenvalues' and    !
    ! 'gamma_eigenvectors' arrays, and will initialize 'n_ions',              !
    ! 'n_branches', 'n_qpts' and 'cell_volume'.                               !
    !=========================================================================!
      implicit none

      ! Dummy variables.
      integer :: phonon_data_unit, ios
      character(len=256) :: line
      integer :: pos, iq, iq_file, im, imode, iv, ia, iatom
      real(kind=dp) :: freq, ir, ref_wgt
      real(kind=dp) :: qa, qb, qc, q_wgt, xr, xi, yr, yi, zr, zi
      logical :: is_gamma, found_gamma

      ! Dummy variables; for testing norm of eigenvectors.
      integer       :: smallest_nrm_dev_iq,smallest_nrm_dev_m,greatest_nrm_dev_iq,greatest_nrm_dev_m
      real(kind=dp) :: tmp,nrm_dev,smallest_nrm_dev,greatest_nrm_dev

      phonon_data_unit = io_file_unit()
      is_gamma = .false.
      found_gamma = .false.

      ! Open seedname.phonon file.
      open(newunit=phonon_data_unit, file=trim(seedname)//".phonon", action='read', iostat=ios)
      if (ios /= 0) CALL io_error("ERROR: Could not open "//trim(seedname)//".phonon")

      ! Parse header; get n_ions, n_branches, n_qpts.
      do
        read(phonon_data_unit, '(A)', iostat=ios) line
        if (ios /= 0) stop "Unexpected EOF in header"

        if (index(line, "Number of ions") > 0) then
          pos = scan(line, "0123456789")
          read(line(pos:), *) n_ions
          ! Allocate these now; will read in atomic posns and species shortly.
          allocate(atomic_positions(3, n_ions))
          allocate(atomic_species(n_ions))
          allocate(atomic_masses(n_ions))

        else if (index(line, "Number of branches") > 0) then
          pos = scan(line, "0123456789")
          read(line(pos:), *) n_branches

        else if (index(line, "Number of wavevectors") > 0) then
          pos = scan(line, "0123456789")
          read(line(pos:), *) n_qpts

        else if (index(line, "Unit cell vectors") > 0) then
          do iv = 1,3
            read(phonon_data_unit, *, iostat=ios) real_lattice(iv,1), real_lattice(iv,2), real_lattice(iv,3)
            if (ios /= 0) CALL io_error("Error reading cell vectors from .phonon file")
          end do

        else if (index(line, "Fractional Co-ordinates") > 0) then
          do ia = 1,n_ions
            read(phonon_data_unit, *, iostat=ios) iatom, atomic_positions(1, ia), atomic_positions(2, ia), &
            atomic_positions(3, ia), atomic_species(ia), atomic_masses(ia)
            if (ios /= 0) CALL io_error("Error reading atomic coordinates from .phonon file")
          end do
        end if

        if (index(line, "END header") > 0) then
          exit
        end if
      end do

      ! Calculate cell volume.
      CALL phonon_eels_calculate_cell_volume(real_lattice,cell_volume)

      ! Calculate RLVs.
      CALL phonon_eels_calculate_RLVs(real_lattice,recip_lattice)

      ! Convert atomic mass from AMU to m_e (natural units).
      do ia=1,n_ions
        atomic_masses(ia) = atomic_masses(ia) * amu_to_me
      end do

      ! Allocate arrays.
      allocate(phonon_eigenvalues(n_branches, n_qpts))
      allocate(phonon_eigenvectors(n_branches, n_ions, 3, n_qpts))
      allocate(qpoint_positions(3, n_qpts))
      allocate(qpoint_weights(n_qpts))
      allocate(gamma_eigenvalues(n_branches))
      allocate(gamma_eigenvectors(n_branches, n_ions, 3))

      ! Loop over q-points, and read-in all phonon data.
      do iq = 1,n_qpts

        ! Read q-point line.
        read(phonon_data_unit, '(A)', iostat=ios) line
        if (ios /= 0) CALL io_error("Error reading q-point line")

        ! Read q-points.
        if (index(line, "q-pt=") == 0) then
          stop "Expected q-pt line"
        end if
        pos = scan(line, "0123456789")
        read(line(pos:), *) iq_file, qa, qb, qc, q_wgt
        qpoint_positions(1, iq) = qa
        qpoint_positions(2, iq) = qb
        qpoint_positions(3, iq) = qc
        qpoint_weights(iq)      = q_wgt

        ! Check if this is Gamma.
        if ( ALL(ABS(qpoint_positions(:,iq)) < 1.0E-6_dp) ) then
          is_gamma = .true.
          found_gamma = .true.
        else
          is_gamma = .false.
        end if

        ! Read frequencies.
        do imode = 1,n_branches
          read(phonon_data_unit, *, iostat=ios) iatom, freq
          if (ios /= 0) stop "Error reading frequencies"
          phonon_eigenvalues(imode, iq) = freq
          if (is_gamma) then
            gamma_eigenvalues(imode) = freq
          end if
        end do

        ! Skip "Phonon Eigenvectors" line.
        read(phonon_data_unit, '(A)') line
        read(phonon_data_unit, '(A)') line

        ! Read eigenvectors.
        ! CASTEP eigenvectors are mass-weighted by default (to see this, test the orthogonality relation).
        ! By dividing by SQRT(mass), we obtain mass-unweighted eigenvectors.
        do im = 1,n_branches
          do ia = 1,n_ions
            read(phonon_data_unit, *, iostat=ios) imode, iatom, xr, xi, yr, yi, zr, zi
            if (ios /= 0) CALL io_error("Error reading eigenvectors from .phonon file")
            phonon_eigenvectors(imode, iatom, 1, iq) = cmplx(xr, xi, dp) / SQRT(atomic_masses(iatom))
            phonon_eigenvectors(imode, iatom, 2, iq) = cmplx(yr, yi, dp) / SQRT(atomic_masses(iatom))
            phonon_eigenvectors(imode, iatom, 3, iq) = cmplx(zr, zi, dp) / SQRT(atomic_masses(iatom))
            if (is_gamma) then
              gamma_eigenvectors(imode, iatom, 1) = cmplx(xr, xi, dp) / SQRT(atomic_masses(iatom))
              gamma_eigenvectors(imode, iatom, 2) = cmplx(yr, yi, dp) / SQRT(atomic_masses(iatom))
              gamma_eigenvectors(imode, iatom, 3) = cmplx(zr, zi, dp) / SQRT(atomic_masses(iatom))
            end if
          end do
        end do
      end do

      close(phonon_data_unit)

      ! Phonon eigenvalues in a .phonon file are in units of inverse cm.
      ! Convert to meV here.
      phonon_eigenvalues = phonon_eigenvalues * inv_cm_to_meV
      gamma_eigenvalues = gamma_eigenvalues * inv_cm_to_meV

      ! Check that Gamma-point data has been found.
      if (.not. found_gamma) then
        CALL io_error("ERROR: The data in the .phonon file does not contain the Gamma-point")
        stop
      end if

      ! For debug: Test normalization of phonon modes (based on iprint level).
      if (iprint > 2) then
        write(stdout,*) ""
        write(stdout,*) "Testing normalization of phonon eigenvectors via:"
        write(stdout,*) ""
        write(stdout,*) "  Norm = sum_{a,i} CONJG(eta^{q,v}_{a,i}) * eta^{q,v}_{a,i} * mass_{a}"
        write(stdout,*) ""
        write(stdout,*) "Where (q,v) indexes the phonon mode, and (a,i) index atom and Cartesian direction"
        write(stdout,*) ""
        write(stdout,*) "For mass-unweighted eigenvectors, this quantity should be exactly 1 for each mode"
        write(stdout,*) ""
        write(stdout,*) "The smallest/largest deviations from unity are shown here:"
        write(stdout,*) ""
        smallest_nrm_dev = 1.00E6_dp
        greatest_nrm_dev = 0.00_dp
        do iq=1,n_qpts
          do imode=1,n_branches
            tmp = 0.0_dp
            do iatom=1,n_ions
              do ia=1,3
                tmp = tmp + CONJG(phonon_eigenvectors(imode,iatom,ia,iq)) &
                * phonon_eigenvectors(imode,iatom,ia,iq) * atomic_masses(iatom)
              end do
            end do
            nrm_dev = ABS(1.00_dp - tmp)
            if ( nrm_dev < smallest_nrm_dev ) then
              smallest_nrm_dev = nrm_dev
              smallest_nrm_dev_iq = iq
              smallest_nrm_dev_m = imode
            end if
            if ( nrm_dev > greatest_nrm_dev ) then
              greatest_nrm_dev = nrm_dev
              greatest_nrm_dev_iq = iq
              greatest_nrm_dev_m = imode
            end if
          end do
        end do
        write(stdout,'(A,I3,A,I5,A,F20.16)') "   Smallest Norm Deviation (mode ",smallest_nrm_dev_m, &
        ", at q-point ",smallest_nrm_dev_iq,") is: ",smallest_nrm_dev
        write(stdout,'(A,I3,A,I5,A,F20.16)') "   Greatest Norm Deviation (mode ",greatest_nrm_dev_m, &
        ", at q-point ",greatest_nrm_dev_iq,") is: ",greatest_nrm_dev
        write(stdout,*) ""
      end if

      ! For debug: Print some important quantities parsed from the .phonon file (based on iprint level).
      if (iprint > 2) then

        write(stdout,*) ""
        write(stdout,*) "The following data was parsed from "//TRIM(seedname)//".phonon:"
        write(stdout,*) ""

        write(stdout,*) "REAL-CELL DATA:"
        do iv=1,3
          write(stdout,*) real_lattice(iv,1),real_lattice(iv,2),real_lattice(iv,3)
        end do
        write(stdout,*) ""

        write(stdout,*) "RECIP-CELL DATA:"
        write(stdout,*) "(calculated from REAL-CELL DATA)"
        do iv=1,3
          write(stdout,*) recip_lattice(iv,1),recip_lattice(iv,2),recip_lattice(iv,3)
        end do
        write(stdout,*) ""

        write(stdout,*) "ATOMIC POSITION DATA:"
        do ia=1,n_ions
          write(stdout,*) "Atom: ",ia,", type: ",atomic_species(ia),",   mass (m_e):",atomic_masses(ia), &
          ", posn = ",atomic_positions(1,ia)," ",atomic_positions(2,ia)," ",atomic_positions(3,ia)
        end do
        write(stdout,*) ""

        write(stdout,*) "QPOINT POSITIONS:"
        do iq=1,n_qpts
          write(stdout,*) "qpt ",iq," = ",qpoint_positions(1,iq)," ",qpoint_positions(2,iq)," ",qpoint_positions(3,iq)
        end do
        write(stdout,*) ""

      end if

      parsed_phonon_file = .true.

    end subroutine phonon_eels_read_phonon_file

    !===========================================================================!
    subroutine phonon_eels_reorder_phonon_bands
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! This subroutine reorders the eigenvectors and eigenvalues in-place,       !
    ! using a band-matching algorithm. This ensures that for quantities         !
    ! that are evaluated along a high-symmetry path, band crossings are         !
    ! correctly detected.                                                       !
    !                                                                           !
    ! Calling this will reorder the 'phonon_eigenvectors' and                   !
    ! 'phonon eigenvalues' arrays in-place.                                     !
    !===========================================================================!
      implicit none

      ! Dummy variables.
      integer                         :: iq, m, n, i, a, best_n
      real(kind=dp)                   :: best_val, dq_nrm
      complex(kind=dp)                :: dot
      logical                         :: just_broke
      real(kind=dp), allocatable      :: overlap(:,:)
      logical, allocatable            :: used(:)
      integer, allocatable            :: perm(:)
      complex(kind=dp), allocatable   :: tmp_vec(:,:,:)
      real(kind=dp), allocatable      :: tmp_val(:)
      real(kind=dp)                   :: dq(3)

      ! Check that eigvecs & eigvals have been parsed already.
      ! If not, then we need to abort because we have no data to reorder.
      if ( (.not. allocated(phonon_eigenvectors)) .or. (.not. allocated(phonon_eigenvalues)) ) then
        CALL io_error("ERROR: Called phonon_eels_reorder_phonon_bands, but eigenvectors and/or eigenvalues have not been parsed")
      end if

      ! Allocate temporary arrays.
      allocate(overlap(n_branches,n_branches))
      allocate(used(n_branches))
      allocate(perm(n_branches))
      allocate(tmp_vec(n_branches,n_ions,3))
      allocate(tmp_val(n_branches))

      ! This is used to detect path breaks.
      just_broke = .false.

      ! Loop over all q-points.
      do iq = 2,n_qpts

          ! Detect path break.
          dq = qpoint_positions(:,iq) - qpoint_positions(:,iq-1)
          dq_nrm = SQRT( dq(1)**2 + dq(2)**2 + dq(3)**2 )
          if (dq_nrm > 0.05) then
            !write(*,*) "Path break at idx: ",iq
            just_broke = .true.
          else
            just_broke = .false.
          end if

          ! Eigenvector overlap computation.
          do m = 1,n_branches
            do n = 1,n_branches
              dot = (0.0_dp, 0.0_dp)
              do a = 1,n_ions
                do i = 1,3
                  dot = dot + CONJG(phonon_eigenvectors(m,a,i,iq-1)) * phonon_eigenvectors(n,a,i,iq)
                end do
              end do
              overlap(m,n) = abs(dot)
            end do
          end do

          ! Only do the matching if there's no path break.
          ! If there is a path break, this is skipped, and the orderings from the previous q-point are used.
          ! This prevents wild discontinuities at path breaks.
          if (.not. just_broke) then
            used = .false.
            do m = 1, n_branches
              best_val = -1.0_dp
              best_n = -1

              do n = 1,n_branches
                if (.not. used(n)) then
                  if (overlap(m,n) > best_val) then
                    best_val = overlap(m,n)
                    best_n = n
                  end if
                end if
              end do

              perm(m) = best_n
              used(best_n) = .true.
            end do
          end if

          ! Reorder into temporary arrays.
          do m = 1,n_branches
            tmp_val(m) = phonon_eigenvalues(perm(m), iq)
            do a = 1,n_ions
              do i = 1,3
                tmp_vec(m,a,i) = phonon_eigenvectors(perm(m),a,i,iq)
              end do
            end do
          end do

          ! Now, swap eigenvalues and eigenvectors in-place.
          phonon_eigenvalues(:,iq) = tmp_val(:)
          phonon_eigenvectors(:,:,:,iq) = tmp_vec(:,:,:)

      end do

      ! Deallocate temporary arrays.
      deallocate(overlap, used, perm, tmp_vec, tmp_val)

    end subroutine phonon_eels_reorder_phonon_bands

    !===========================================================================!
    subroutine phonon_eels_detect_hsps
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                         !
    !                                                                           !
    ! If the phonons are supplied along a dispersion path, this subroutine      !
    ! will attempt to detect the high-symmetry points along that path, by       !
    ! tracking changes in the path gradient and/or locating path breaks. This   !
    ! is useful, because CASTEP does not write the location of the              !
    ! high-symmetry points to the .phonon file.                                 !
    !                                                                           !
    ! The locations (q-point indexes) of the detected high-symmetry points,     !
    ! along with the q-point spacings (in inverse Angstrom) between those       !
    ! points, are then written in the header of any output files that contain   !
    ! data evaluated along the dispersion (impact_eels_intensity, affs, dwfs,   !
    ! etc). This is helpful for plotting purposes.                              !
    !                                                                           !
    ! Calling this will allocate & populate the 'hsp_idxs' and                  !
    ! 'path_spacing_norms' arrays.                                              !
    !===========================================================================!
      implicit none

      ! Dummy variables.
      integer                         :: iq,i,N_unique_hsps
      real(kind=dp)                   :: dq(3), dq_old(3), qx, qy, qz, path_spacing_norm
      integer                         :: hsp_idxs_buffer(50),unique_hsp_idxs_buffer(50)

      ! Check that eigvecs & eigvals have been parsed already.
      ! If not, then we need to abort because we have no data from which to detect HSPs.
      if ( (.not. allocated(phonon_eigenvectors)) .or. (.not. allocated(phonon_eigenvalues)) ) then
        CALL io_error("ERROR: Called phonon_eels_detect_hsps, but eigenvectors and/or eigenvalues have not been parsed")
      end if

      ! We always start the path with a HSP.
      N_hsps = 1
      hsp_idxs_buffer(N_hsps) = 1

      ! Now, loop over all q-points, and detect HSPs.
      dq_old = qpoint_positions(:,2) - qpoint_positions(:,1)
      do iq=3,n_qpts
        dq = qpoint_positions(:,iq) - qpoint_positions(:,iq-1)
        ! Case of a regular HSP location, where the path gradient changes discontinuously.
        if ( ANY(ABS(dq - dq_old) > 1.0E-5_dp) ) then
          N_hsps = N_hsps + 1
          hsp_idxs_buffer(N_hsps) = iq - 1
          dq_old = dq
          cycle
        end if
        ! Special case of a path break where the gradient happens to be the same on both sides.
        ! (the above check would otherwise miss this)
        if ( ANY(ABS(dq) > 1.0E-1_dp) ) then
          N_hsps = N_hsps + 1
          hsp_idxs_buffer(N_hsps) = iq - 1
          dq_old = dq
          cycle
        end if
        dq_old = dq
      end do

      ! And we always finish the path with a HSP.
      N_hsps = N_hsps + 1
      hsp_idxs_buffer(N_hsps) = n_qpts

      ! Now, count unique HSPs.
      ! (the above code will count a path break as two HSPs, so we need to adjust this if necessary)
      N_unique_hsps = 1
      unique_hsp_idxs_buffer(1) = hsp_idxs_buffer(1)
      do i=2,N_hsps
        if ( (hsp_idxs_buffer(i) - hsp_idxs_buffer(i-1)) .eq. 1 ) then
          cycle
        else
          N_unique_hsps = N_unique_hsps + 1
          unique_hsp_idxs_buffer(N_unique_hsps) = hsp_idxs_buffer(i)
        end if
      end do

      ! Allocate 'hsp_idxs' and store the indexes in there.
      allocate(hsp_idxs(N_unique_hsps))
      do i=1,N_unique_hsps
        hsp_idxs(i) = unique_hsp_idxs_buffer(i)
      end do

      ! Update N_hsps.
      N_hsps = N_unique_hsps

      ! Allocate path_spacing_norms.
      allocate(path_spacing_norms(N_hsps-1))

      ! Get spacing of q-points between the HSPs.
      do i=1,N_hsps-1
        iq = hsp_idxs(i)
        dq = qpoint_positions(:,iq+2) - qpoint_positions(:,iq+1)
        ! Get the magnitude of this q-vector, in units of inverse Angstrom.
        qx = dq(1) * recip_lattice(1,1) + dq(2) * recip_lattice(2,1) + dq(3) * recip_lattice(3,1)
        qy = dq(1) * recip_lattice(1,2) + dq(2) * recip_lattice(2,2) + dq(3) * recip_lattice(3,2)
        qz = dq(1) * recip_lattice(1,3) + dq(2) * recip_lattice(2,3) + dq(3) * recip_lattice(3,3)
        path_spacing_norm = SQRT( qx**2 + qy**2 + qz**2 )
        path_spacing_norms(i) = path_spacing_norm
      end do

      detected_hsps = .true.

    end subroutine phonon_eels_detect_hsps

    !===========================================================================!
    subroutine phonon_eels_write_hsps_header(unit_number)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, May 2026                                         !
    !                                                                           !
    ! This subroutine writes a header containing the locations (q-point idxs)   !
    ! of the HSPs to the supplied unit number, along with the q-point spacings  !
    ! (in inverse Angstrom) between those HSPs. It can be called by any         !
    ! subroutine that requires HSP data to be written in the output file        !
    ! header.                                                                   !
    !===========================================================================!
      implicit none

      ! Subroutine arguments.
      integer,INTENT(IN) :: unit_number

      ! Dummy variables.
      integer :: i

      ! Check that we actually have the HSPs first.
      if (.not. detected_hsps) then
        CALL io_error("ERROR: Cannot write high-symmetry points header, because HSPs have not been detected")
      end if

      ! Write the HSPs header to the supplied unit_number.
      write(unit_number,'(A)') "  High-symmetry points and/or path breaks detected at the following q_idxs:"
      do i=1,N_hsps
        write(unit_number,'(I5,A)',advance="no") hsp_idxs(i)," "
      end do
      write(unit_number,'(A)') ""
      write(unit_number,'(A)') ""
      write(unit_number,'(A)') "  The path spacings (in inverse Angstrom) between these points are:"
      do i=1,N_hsps-1
        write(unit_number,'(A,ES10.4)',advance="no") "    ",path_spacing_norms(i)
      end do
      write(unit_number,'(A)') ""

    end subroutine phonon_eels_write_hsps_header

    !===========================================================================!
    subroutine phonon_eels_calculate_RLVs(L,L_rlv)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! Given a real-space lattice matrix 'L', this subroutine calculates the     !
    ! RLVs, and stores them in 'L_rlv'.                                         !
    !===========================================================================!
      implicit none

      ! Subroutine arguments.
      real(kind=dp),INTENT(IN)  :: L(3,3)
      real(kind=dp),INTENT(OUT) :: L_rlv(3,3)

      ! Dummy variables.
      real(kind=dp)              :: L_inv(3,3), V

      ! Get cell volume (i.e. determinant of L).
      CALL phonon_eels_calculate_cell_volume(L,V)

      if (ABS(V) < 1.0E-12_dp) then
        CALL io_error("ERROR in routine phonon_eels_calculate_RLVs: Lattice matrix is singular")
      end if

      L_inv(1,1) =  (L(2,2)*L(3,3) - L(2,3)*L(3,2)) / V
      L_inv(1,2) = -(L(1,2)*L(3,3) - L(1,3)*L(3,2)) / V
      L_inv(1,3) =  (L(1,2)*L(2,3) - L(1,3)*L(2,2)) / V

      L_inv(2,1) = -(L(2,1)*L(3,3) - L(2,3)*L(3,1)) / V
      L_inv(2,2) =  (L(1,1)*L(3,3) - L(1,3)*L(3,1)) / V
      L_inv(2,3) = -(L(1,1)*L(2,3) - L(1,3)*L(2,1)) / V

      L_inv(3,1) =  (L(2,1)*L(3,2) - L(2,2)*L(3,1)) / V
      L_inv(3,2) = -(L(1,1)*L(3,2) - L(1,2)*L(3,1)) / V
      L_inv(3,3) =  (L(1,1)*L(2,2) - L(1,2)*L(2,1)) / V

      L_rlv = twopi * TRANSPOSE(L_inv)

    end subroutine phonon_eels_calculate_RLVs

    !===========================================================================!
    subroutine phonon_eels_calculate_cell_volume(L,V)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, April 2026                                       !
    !                                                                           !
    ! Given a real-space lattice matrix 'L', this subroutine calculates the     !
    ! volume of the unit cell, and stores it in 'V'.                            !
    !===========================================================================!
      implicit none

      ! Subroutine arguments.
      real(kind=dp),INTENT(IN)  :: L(3,3)
      real(kind=dp),INTENT(OUT) :: V

      V = L(1,1)*(L(2,2)*L(3,3) - L(2,3)*L(3,2)) - L(1,2)*(L(2,1)*L(3,3) &
          - L(2,3)*L(3,1)) + L(1,3)*(L(2,1)*L(3,2) - L(2,2)*L(3,1))

    end subroutine phonon_eels_calculate_cell_volume

    !===========================================================================!
    subroutine phonon_eels_get_normal_rotation_matrix(surface_plane,R)
    !===========================================================================!
    ! J. A. J. Whaley-Baldwin, June 2026                                        !
    !                                                                           !
    ! This subroutine constructs an orthogonal matrix that aligns the chosen    !
    ! reciprocal lattice axis ('normal_axis') with the z-direction, and         !
    ! returns the resulting matrix in 'R'.                                      !
    !                                                                           !
    ! This is useful for the aloof theory, since the z-direction is always      !
    ! defined as the normal to the crystal surface. So, for non-cubic crystals, !
    ! different surface normals give a different result.                        !
    !                                                                           !
    ! For cases where the crystal is already oriented along the chosen surface  !
    ! normal, this subroutine just returns the identity matrix.                 !
    !                                                                           !
    ! Once constructed, any vector/tensor quantities can then be rotated like:  !
    !                                                                           !
    !        v_rot       =  matmul(R,v)                                         !
    !        tensor_rot  =  matmul( R, matmul(tensor_rot,TRANSPOSE(R)) )        !
    !                                                                           !
    !===========================================================================!
      implicit none

      ! Subroutine arguments.
      character(len=*),INTENT(IN)  :: surface_plane
      real(kind=dp),INTENT(OUT)    :: R(3,3)

      ! Dummy variables.
      real(dp) :: u(3), z(3)
      real(dp) :: v(3), K(3,3)
      real(dp) :: s, c

      ! We need the cell data from the .phonon file first.
      if (.not. parsed_phonon_file) then
        CALL phonon_eels_read_phonon_file
      end if

      ! Get the appropriate RLV.
      ! The rotation matrix will then be constructed so that this RLV lies along z.
      if (index(surface_plane,"ab") > 0) then
        u = recip_lattice(3,:)
      else if (index(surface_plane,"ac") > 0) then
        u = recip_lattice(2,:)
      else if (index(surface_plane,"bc") > 0) then
        u = recip_lattice(1,:)
      else
        CALL io_error("ERROR: Crystal surface plane for aloof vib-EELS must be one of 'ab', 'ac' or 'bc'")
      end if

      ! Input vector.
      u = u / norm2(u)

      ! Target direction for the preferential axis.
      z = [0.0_dp, 0.0_dp, 1.0_dp]

      ! If u is already aligned with z, then just return identity.
      if ( ALL( ABS(u-z) < 1.0E-6_dp) ) then
        R = 0.0_dp
        R(1,1) = 1.0_dp
        R(2,2) = 1.0_dp
        R(3,3) = 1.0_dp
        return
      end if

      ! Rotation axis.
      v(1) = u(2)*z(3) - u(3)*z(2)
      v(2) = u(3)*z(1) - u(1)*z(3)
      v(3) = u(1)*z(2) - u(2)*z(1)

      s = norm2(v)
      c = dot_product(u,z)

      ! Skew-symmetric matrix.
      K = 0.0_dp
      K(1,2) = -v(3)
      K(1,3) =  v(2)
      K(2,1) =  v(3)
      K(2,3) = -v(1)
      K(3,1) = -v(2)
      K(3,2) =  v(1)

      ! Identity.
      R = 0.0_dp
      R(1,1) = 1.0_dp
      R(2,2) = 1.0_dp
      R(3,3) = 1.0_dp

      ! Rodrigues formula.
      R = R + K + matmul(K,K)*(1.0_dp-c) / (s*s)

    end subroutine phonon_eels_get_normal_rotation_matrix

end module od_phonon_eels