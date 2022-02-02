program bands
  implicit none

  integer, parameter:: dp = selected_real_kind(15, 300)  ! Define double precision
  integer, parameter:: num_kpoints = 1                ! Number of kpoints
  integer, parameter:: num_spins = 2                  ! Number of spins
  integer, parameter:: max_eigenv = 1             ! Number of bands included in matrix elements
  integer, parameter:: num_electrons(num_spins) = 1   ! Number of electrons of each spin
  integer, parameter:: num_eigenvalues(num_spins) = 1 ! Number of eigenvalues of each spin
  real(dp):: efermi(num_spins) ! Fermi level calculated from electronic
  ! structure code for each spin
  real(dp):: cell_lattice(1:3, 1:3) ! Reals-space lattice vectors
  real(dp):: kpoint_positions(num_kpoints, 1:3) ! K-point position vectors in
  ! fractions of Brillouin Zone
  real(dp):: kpoint_weight(num_kpoints)        ! K-point weight
  real(dp):: band_energy(max_eigenv, num_spins, num_kpoints) ! Energy
  ! Eigenvalue matrix

  integer :: ik, is, ib, i, jb, band_unit = 6, orb, indx
  open (unit=band_unit, form='formatted', file="band.out")

  write (band_unit, '(a,i4)') 'Number of k-points ', num_kpoints
  write (band_unit, '(a,i1)') 'Number of spin components ', num_spins
  if (num_spins == 2) then
    write (band_unit, '(a,2g10.4)'), 'Number of electrons ', num_electrons(:)
  else
    write (band_unit, '(a,g10.4)'), 'Number of electrons ', num_electrons(:)
  end if
  if (num_spins == 2) then
    write (band_unit, '(a,2i4)'), 'Number of eigenvalues ', num_eigenvalues(:)
  else
    write (band_unit, '(a,i4)'), 'Number of eigenvalues ', num_eigenvalues(:)
  end if
  if (num_spins == 2) then
    write (band_unit, '(a,2f12.6)') 'Fermi energy (in atomic units) ', efermi(:)
  else
    write (band_unit, '(a,f12.6)') 'Fermi energies (in atomic units) ', efermi(:)
  end if
  write (band_unit, '(a)') 'Unit cell vectors'
  write (band_unit, '(3f12.6)') cell_lattice(1, :)
  write (band_unit, '(3f12.6)') cell_lattice(2, :)
  write (band_unit, '(3f12.6)') cell_lattice(3, :)

  do ik = 1, num_kpoints
    write (band_unit, '(a,i4,4f12.8)') 'K-point ', ik, kpoint_positions(num_kpoints, 1:3), &
           & kpoint_weight(num_kpoints)
    do is = 1, num_spins
      write (band_unit, '(a,i1)') 'Spin component ', is
      do ib = 1, num_eigenvalues(is)
        write (band_unit, *) band_energy(ib, is, ik)
      end do
    end do
  end do

end program bands
