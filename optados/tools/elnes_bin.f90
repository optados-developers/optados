program elnes_bin
  implicit none

  integer, parameter:: dp = selected_real_kind(15, 300)
  integer, parameter:: tot_core_projectors = 1  ! Total number of core states included in matrix elements
  integer, parameter:: max_eigenvalues = 1      ! Number of bands included in matrix elements
  integer, parameter:: num_kpoints = 1          ! Number of k-points
  integer, parameter:: num_spins = 1            ! Number of spins
  integer, parameter:: num_eigenvalues(1:num_spins) = 1    ! Number of eigenvalues
  integer:: species(1:tot_core_projectors) ! Atomic species associated with each projector
  integer:: ion(1:tot_core_projectors)     ! Ion associated with each projector
  integer:: n(1:tot_core_projectors)       ! Principal quantum number associated with each projector
  integer:: lm(1:tot_core_projectors)      ! Angular momentum quantum number associated with each projector
  real(dp):: elnes_mat(tot_core_projectors, max_eigenvalues, 3, num_kpoints, num_spins) ! matrix elements

  integer :: nk, ns, nb, i, jb, elnes_file = 6, orb, indx
  open (unit=elnes_file, form='unformatted', file="elnes.out")

  write (elnes_file) tot_core_projectors
  write (elnes_file) max_eigenvalues
  write (elnes_file) num_kpoints
  write (elnes_file) num_spins
  write (elnes_file) species(1:tot_core_projectors)
  write (elnes_file) ion(1:tot_core_projectors)
  write (elnes_file) n(1:tot_core_projectors)
  write (elnes_file) lm(1:tot_core_projectors)

  do nk = 1, num_kpoints
    do ns = 1, num_spins
      write (elnes_file) (((elnes_mat(orb, nb, indx, nk, ns), orb=1,&
           &tot_core_projectors), nb=1, num_eigenvalues(ns)), indx=1, 3)
    end do
  end do

end program elnes_bin
