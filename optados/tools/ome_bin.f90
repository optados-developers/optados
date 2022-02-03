program ome_bin
  implicit none
  integer, parameter:: dp = selected_real_kind(15, 300) ! Define double precision
  real(dp):: file_version = 1.0_dp          ! File version
  character(len=80):: file_header         ! File header comment
  integer, parameter:: num_kpoints = 1                   ! Number of k-points
  integer, parameter:: num_spins = 1                     ! Number of spins
  integer, parameter:: max_eigenv = 1             ! Number of bands included in matrix elements
  integer :: num_eigenvalues(1:num_spins)   ! Number of eigenvalues per spin channel
  complex(dp):: optical_mat(max_eigenv, max_eigenv, 1:3, 1:num_kpoints, num_spins) ! OMEs

  integer :: ik, is, ib, i, jb, ome_unit = 6

  open (unit=ome_unit, form='unformatted', file="ome.out")

  write (ome_unit) file_version
  write (ome_unit) file_header

  do ik = 1, num_kpoints
    do is = 1, num_spins
      write (ome_unit) (((optical_mat(ib, jb, i, ik, is), ib=1, num_eigenvalues(is)), &
           &jb=1, num_eigenvalues(is)), i=1, 3)
    end do
  end do

end program ome_bin
