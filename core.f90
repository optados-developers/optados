!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 


module od_core

  use od_constants, only : dp
  implicit none
  private
  public :: core_calculate

  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_dos



contains

  subroutine core_calculate
    use od_electronic, only  : elec_read_elnes_mat, elnes_mat,  elnes_mwab
    use od_dos_utils, only : dos_utils_calculate
    use od_comms, only : on_root

    implicit none

    ! read in the core matrix elements from disk
    call  elec_read_elnes_mat
!    (elnes_mat(orb,nb,indx,nk,ns),indx=1,3)

    call core_prepare_matrix_elements

    call dos_utils_calculate(matrix_weights, weighted_dos)

    if (on_root) then
       call write_core
    endif

  end subroutine core_calculate

  ! Private routines

  subroutine core_prepare_matrix_elements
    use od_electronic, only  : elnes_mat,  elnes_mwab, nbands, nspins,num_electrons, electrons_per_state
    use od_comms, only : on_root, my_node_id
    use od_cell, only : nkpoints, cell_volume, num_kpoints_on_node
    use od_parameters, only : optics_geom, optics_qdir
    use od_io, only : io_error

    real(kind=dp), dimension(3) :: qdir 
    real(kind=dp) :: q_weight
    integer :: N, N_spin, n_eigen, orb
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp) :: g

    num_occ = 0.0_dp
    do N_spin=1,nspins
       num_occ(N_spin) = num_electrons(N_spin)
    enddo

    if (electrons_per_state==2) then 
       num_occ(1) = num_occ(1)/2.0_dp
    endif

    allocate(matrix_weights(elnes_mwab%norbitals,elnes_mwab%nbands,num_kpoints_on_node(my_node_id),nspins))
    matrix_weights=0.0_dp

    qdir=optics_qdir
    q_weight=((qdir(1)**2.0_dp)+(qdir(2)**2.0_dp)+(qdir(3)**2.0_dp))**0.5_dp
    if(q_weight<0.001_dp)&
         call io_error("Error:  please check optics_qdir, norm close to zero")


    do N=1,num_kpoints_on_node(my_node_id)                      ! Loop over kpoints
       do N_spin=1,nspins                                    ! Loop over spins
             do n_eigen=(nint(num_occ(N_spin))+1),nbands    ! Loop over unoccupied states
                do orb=1,elnes_mwab%norbitals
!                g = (((qdir(1)*elnes_mat(orb,n_eigen,1,N,N_spin))+ &
!                     (qdir(2)*elnes_mat(orb,n_eigen,2,N,N_spin))+ &
!                     (qdir(3)*elnes_mat(orb,n_eigen,3,N,N_spin)))/q_weight)
                   matrix_weights(orb,n_eigen,N,N_spin) = ( &
                        elnes_mat(orb,n_eigen,1,N,N_spin)*conjg(elnes_mat(orb,n_eigen,1,N,N_spin)) + &
                        elnes_mat(orb,n_eigen,2,N,N_spin)*conjg(elnes_mat(orb,n_eigen,2,N,N_spin)) + &
                        elnes_mat(orb,n_eigen,3,N,N_spin)*conjg(elnes_mat(orb,n_eigen,3,N,N_spin)) ) /3.0_dp
!                matrix_weights(orb,n_eigen,N,N_spin) = real(g*conjg(g),dp)
                !           matrix_weights(n_eigen,n_eigen2,N,N_spin) = 1.0_dp  ! 
             end do
          end do
       end do
    end do

  end subroutine core_prepare_matrix_elements

  subroutine write_core
    !***************************************************************
    ! This subroutine writes out the Core loss function

    use od_constants, only : bohr2ang
    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : dos_nbins
    use od_electronic, only: nbands, num_electrons, nspins, elnes_mwab
    use od_io, only : seedname, io_file_unit
    use od_dos_utils, only : E

    integer :: N
    real(kind=dp) ::dE
    integer :: core_unit,orb

    dE=E(2)-E(1)

    ! Open the output file
    core_unit = io_file_unit()
    open(unit=core_unit,action='write',file=trim(seedname)//'.core')

    ! Write into the output file
    write(core_unit,*)'#*********************************************'
    write(core_unit,*)'#            Core loss function               '
    write(core_unit,*)'#*********************************************'

    weighted_dos=weighted_dos*bohr2ang**2


    do orb=1,elnes_mwab%norbitals
       DO N=1,dos_nbins
          write(core_unit,*)E(N),weighted_dos(N,1,orb)
       END DO
    end do

    close(core_unit)


  end subroutine write_core

end module od_core
