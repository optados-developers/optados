program pdos_bin
implicit none
integer,parameter:: dp=selected_real_kind(15,300) ! Define double precision
integer,parameter:: num_kpoints=1           ! Number of k-points
integer,parameter:: num_spins=1             ! Number of spins
integer,parameter:: num_popn_orb=1          ! Number of LCAO projectors
integer,parameter:: max_eigenv=1       ! Number of bands included in matrix elements
real(dp):: file_version=1.0_dp  ! File version
integer:: species(1:num_popn_orb)! Atomic species associated with each projector
integer:: ion(1:num_popn_orb)    ! Ion associated with each projector
integer:: am_channel(1:num_popn_orb)     ! Angular momentum channel
integer:: num_eigenvalues(1:num_spins)   ! Number of eigenvalues per spin channel
real(dp):: kpoint_positions(1:num_kpoints,1:3) ! k_x, k_y, k_z in fractions of BZ
real(dp):: pdos_weights(1:num_popn_orb,max_eigenv,num_kpoints,num_spins)!Matrix elements
character(len=80):: file_header ! File header comment


integer :: nk,ns,nb, pdos_file=6

open(unit=pdos_file, form='unformatted', file="pdos.out")

write(pdos_file) file_header
write(pdos_file) file_version
write(pdos_file) num_kpoints
write(pdos_file) num_spins
write(pdos_file) num_popn_orb
write(pdos_file) max_eigenv
write(pdos_file) species(1:num_popn_orb)
write(pdos_file) ion(1:num_popn_orb)
write(pdos_file) am_channel(1:num_popn_orb)
       
do nk=1,num_kpoints
 write(pdos_file) nk, kpoint_positions(nk,:)
 do ns = 1,num_spins
  write(pdos_file) ns
  write(pdos_file) num_eigenvalues(ns)
  do nb = 1,num_eigenvalues(ns)
   write(pdos_file) real(pdos_weights(num_popn_orb,nb,nk,ns))
  end do
 end do
end do


end program pdos_bin
