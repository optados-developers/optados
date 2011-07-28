!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
module od_optics

  implicit none
  private
  public :: optics_calculate

  type :: graph_labels
     character(20) :: name
     character(40) :: title
     character(20) :: x_label
     character(20) :: y_label
     character(20) :: legend_a
     character(20) :: legend_b
  end type graph_labels


  integer, parameter :: dp=selected_real_kind(15,300)

  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: dos_matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_jdos
  real(kind=dp),allocatable, public, dimension(:,:) :: weighted_dos_at_e
  real(kind=dp),allocatable, public, dimension(:,:) :: dos_at_e

  real(kind=dp),allocatable, dimension(:,:,:,:) :: epsilon
  real(kind=dp),allocatable, dimension(:,:) :: conduct
  real(kind=dp),allocatable, dimension(:,:) :: refract
  real(kind=dp),allocatable, dimension(:,:) :: loss_fn
  real(kind=dp),allocatable, dimension(:) :: absorp
  real(kind=dp),allocatable, dimension(:) :: reflect

  real(kind=dp),allocatable, dimension(:) :: intra
  real(kind=dp) :: q_weight 
  real(kind=dp) :: N_eff
  real(kind=dp) :: N_eff2
  real(kind=dp) :: N_eff3
  integer :: N_geom
  real(kind=dp) :: e_fermi
  integer :: N
  integer :: N2


  real(kind=dp), parameter :: epsilon_0 = 8.8541878176E-12  ! need to put in correct number
  real(kind=dp), parameter :: e_charge =  1.60217646E-19 ! need to put in correct number 
  real(kind=dp), parameter :: e_mass =  9.10938E-31 ! need to put in correct number
  real(kind=dp), parameter :: hbar =  1.054571628E-34 ! need to put in correct number 
  real(kind=dp), parameter :: c_speed =  299792458 ! need to put in correct number

contains

  subroutine optics_calculate
    !
    !  Program to calculate optical properties
    !

    use od_constants, only : dp
    use od_electronic, only : band_gradient, elec_read_band_gradient, nbands, nspins
    use od_cell, only : cell_volume, num_kpoints_on_node 
    use od_jdos_utils, only : jdos_utils_calculate
    use od_comms, only : on_root, my_node_id
    use od_parameters, only : optics_geom, compute_efermi, adaptive, linear, fixed, optics_intraband, &
         optics_drude_broadening 
    use od_dos_utils, only : dos_utils_calculate_at_e, efermi_fixed, efermi_adaptive, efermi_linear
    use od_io, only : stdout 

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+=============================== Optics Calculation =========================+'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif

    ! Get Ef
    if(compute_efermi) then !this is assuming that I have already run the DOS!!!!  
       if(adaptive) e_fermi = efermi_adaptive
       if(linear) e_fermi = efermi_linear
       if(fixed)  e_fermi = efermi_fixed
       !   else call io_error ("OPTICS: No Ef set")
    endif

    ! Get information from .cst_ome file 
    call elec_read_band_gradient

    ! Form matrix element
    call make_weights

    ! Send matrix element to jDOS routine and get weighted jDOS back
    call jdos_utils_calculate(matrix_weights, weighted_jdos)

    ! Calculate weighted DOS at Ef for intraband term
    if(optics_intraband)then 
       allocate(dos_matrix_weights(size(matrix_weights,5),nbands,num_kpoints_on_node(my_node_id),nspins))
       allocate(dos_at_e(3,nspins))
       allocate(weighted_dos_at_e(nspins,size(matrix_weights,5)))
       weighted_dos_at_e = 0.0_dp
       do N=1,size(matrix_weights,5)
          do N2=1,nbands    
             dos_matrix_weights(N,N2,:,:) = matrix_weights(N2,N2,:,:,N)  
          enddo
       enddo
       call dos_utils_calculate_at_e(e_fermi,dos_matrix_weights,weighted_dos_at_e,dos_at_e) 
    endif

    if(on_root) then
       ! Calculate epsilon_2
       call calc_epsilon_2

       ! Calculate epsilon_1
       call calc_epsilon_1

       ! Calculate other optical properties
       if (.not. index(optics_geom,'tensor')>0)then
          call calc_conduct
          call calc_refract
          call calc_loss_fn
          call calc_absorp
          call calc_reflect
       end if

       ! Write everything out
       call write_epsilon
       if (.not. index(optics_geom,'tensor')>0)then
          call write_conduct
          call write_refract
          call write_loss_fn
          call write_absorp
          call write_reflect
       end if
    endif

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+============================== Optics Calculation End ======================+'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif

  end subroutine optics_calculate

  ! Subroutines go here 

  !***************************************************************
  subroutine make_weights
    !***************************************************************
    use od_constants, only : dp
    use od_electronic, only : nbands, nspins, band_gradient, num_electrons, &
         electrons_per_state, band_energy
    use od_cell, only : nkpoints, cell_volume, num_kpoints_on_node, cell_get_symmetry, &
         num_crystal_symmetry_operations, crystal_symmetry_operations
    use od_parameters, only : optics_geom, optics_qdir, legacy_file_format, scissor_op
    use od_io, only : io_error
    use od_comms, only : my_node_id

    real(kind=dp), dimension(3) :: qdir 
    real(kind=dp), dimension(3) :: qdir1 
    real(kind=dp), dimension(3) :: qdir2
    real(kind=dp) :: q_weight1 
    real(kind=dp) :: q_weight2 
    integer :: N, i, j
    integer :: N_in
    integer :: N_spin
    integer :: N2, N3
    integer :: n_eigen
    integer :: n_eigen2
    integer :: num_symm
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp), dimension(3) :: g
    real(kind=dp) :: factor

    num_symm=0
    if (.not.legacy_file_format) then 
       call cell_get_symmetry
       num_symm = num_crystal_symmetry_operations

       !      Do N=1,num_symm
       !         Do N2=1,3
       !            print *, crystal_symmetry_operations(1,N2,N), crystal_symmetry_operations(2,N2,N), crystal_symmetry_operations(3,N2,N)
       !         end do
       !         print *,'# '
       !      end do
    end if

    num_occ = 0.0_dp
    do N_spin=1,nspins
       num_occ(N_spin) = num_electrons(N_spin)
    enddo

    if (electrons_per_state==2) then 
       num_occ(1) = num_occ(1)/2.0_dp
    endif

    if (.not. index(optics_geom,'tensor')>0) then ! I can rewrite this in a simplier way??
       N_geom=1
    elseif (index(optics_geom,'tensor')>0) then 
       N_geom=6
    endif

    allocate(matrix_weights(nbands,nbands,num_kpoints_on_node(my_node_id),nspins,N_geom))
    matrix_weights=0.0_dp

    if (index(optics_geom,'polar')>0) then 
       qdir=optics_qdir
       q_weight=((qdir(1)**2)+(qdir(2)**2)+(qdir(3)**2))**0.5_dp
       if(q_weight<0.001_dp)&
            call io_error("Error:  please check optics_qdir, norm close to zero")
    endif

    if (index(optics_geom,'unpolar')>0) then 
       if(optics_qdir(3)==0)then 
          qdir1(1)=0.0_dp
          qdir1(2)=0.0_dp
          qdir1(3)=1.0_dp
       else
          qdir1(1)=1.0_dp
          qdir1(2)=1.0_dp
          qdir1(3)=-(optics_qdir(1)+optics_qdir(2))/optics_qdir(3)
       endif
       qdir2(1)=(optics_qdir(2)*qdir1(3))-(optics_qdir(3)*qdir1(2))            
       qdir2(2)=(optics_qdir(3)*qdir1(1))-(optics_qdir(1)*qdir1(3))            
       qdir2(3)=(optics_qdir(1)*qdir1(2))-(optics_qdir(2)*qdir1(1))  
       q_weight1=((qdir1(1)**2)+(qdir1(2)**2)+(qdir1(3)**2))**0.5_dp
       q_weight2=((qdir2(1)**2)+(qdir2(2)**2)+(qdir2(3)**2))**0.5_dp
    endif

    N_in = 1  ! 0 = no inversion, 1 = inversion   
    g = 0.0_dp

    do N=1,num_kpoints_on_node(my_node_id)                   ! Loop over kpoints
       do N_spin=1,nspins                                    ! Loop over spins
          do n_eigen=1,nbands                                ! Loop over state 1 
             do n_eigen2=n_eigen,nbands                      ! Loop over state 2  
                if(band_energy(n_eigen,N_spin,N)>e_fermi.and.n_eigen/=n_eigen2) cycle
                if(band_energy(n_eigen2,N_spin,N)<e_fermi.and.n_eigen/=n_eigen2) cycle
                factor = 0.0_dp
                if(n_eigen2==n_eigen)then 
                   factor = 1.0_dp
                elseif(band_energy(n_eigen2,N_spin,N)>e_fermi)then 
                   factor = 1.0_dp/((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)&
                        +scissor_op)**2)  
                endif
                if(index(optics_geom,'unpolar')>0)then
                   if (num_symm==0) then 
                      g(1) = (((qdir1(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir1(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir1(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight1)
                      g(2) = (((qdir2(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir2(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir2(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = &
                           0.5_dp*factor*(real(g(1)*conjg(g(1)),dp)+real(g(2)*conjg(g(2)),dp))
                   else
                      do N2=1,num_symm
                         do N3=1,1+N_in
                            do i=1,3
                               qdir(i)=0.0_dp
                               do j=1,3
                                  qdir(i) = qdir(i) + ((-1.0_dp)**(N3+1))*&
                                       (crystal_symmetry_operations(j,i,N2)*qdir1(j))
                               end do
                            end do
                            g(1)=(((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight1)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) =  &
                                 matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) + &
                                 (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)*factor
                            g(1) = 0.0_dp
                            do i=1,3 ! if I include an extra variable I can merge this and the last do loops
                               qdir(i)=0.0_dp
                               do j=1,3
                                  qdir(i) = qdir(i) + ((-1.0_dp)**(N3+1))*&
                                       (crystal_symmetry_operations(j,i,N2)*qdir2(j))
                               end do
                            end do
                            g(1)=(((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight2)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) =  &
                                 matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) + &
                                 (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)*factor
                         end do
                      end do
                   end if
                elseif(index(optics_geom,'polar')>0)then
                   if (num_symm==0) then 
                      g(1) = (((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = factor*real(g(1)*conjg(g(1)),dp)
                      !           matrix_weights(n_eigen,n_eigen2,N,N_spin) = 1.0_dp  ! for JDOS 
                   else
                      do N2=1,num_symm
                         do N3=1,1+N_in
                            do i=1,3
                               qdir(i)=0.0_dp
                               do j=1,3
                                  qdir(i) = qdir(i) + ((-1.0_dp)**(N3+1))*&
                                       (crystal_symmetry_operations(j,i,N2)*optics_qdir(j))
                               end do
                            end do
                            g(1)=0.0_dp
                            g(1)=(((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) =  &
                                 matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) + &
                                 (1.0_dp/Real((num_symm*(N_in+1)),dp))*factor*real(g(1)*conjg(g(1)),dp)
                         end do
                      end do
                   end if
                elseif(index(optics_geom,'poly')>0)then
                   if (num_symm==0) then 
                      do N2=1,3  
                         g(N2) = band_gradient(n_eigen,n_eigen2,N2,N,N_spin)
                         !if (n_eigen==1) then 
                         !if (n_eigen2==NINT(num_occ(N_spin)+1)) then 
                         !print *, N, N2 
                         !print *, g(N2), abs(g(N2))
                         !print *, band_gradient(n_eigen2,n_eigen,N2,N,N_spin)
                         !print *, band_gradient(n_eigen,n_eigen,N2,N,N_spin), band_gradient(n_eigen2,n_eigen2,N2,N,N_spin)
                         !end if 
                         !end if 
                      end do
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = (factor/3.0_dp)*&
                           (real(g(1)*conjg(g(1)),dp)+ real(g(2)*conjg(g(2)),dp) + real(g(3)*conjg(g(3)),dp))
                      !                 print *, n_eigen, n_eigen2, N, matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom)
                      !                 print *, band_energy(n_eigen2,N_spin,N), band_energy(n_eigen,N_spin,N)
                   else
                      do N2=1,num_symm
                         do N3=1,1+N_in
                            qdir = 0.0_dp
                            qdir1=0.0_dp
                            qdir2=0.0_dp
                            do i=1,3
                               qdir(i) = ((-1.0_dp)**(N3+1))*crystal_symmetry_operations(1,i,N2)
                               qdir1(i) = ((-1.0_dp)**(N3+1))*crystal_symmetry_operations(2,i,N2)
                               qdir2(i) = ((-1.0_dp)**(N3+1))*crystal_symmetry_operations(3,i,N2)
                            end do
                            g = 0.0_dp
                            g(1)=((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))
                            g(2)=((qdir1(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir1(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir1(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))
                            g(3)=((qdir2(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                                 (qdir2(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                                 (qdir2(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) =  &
                                 matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) + &
                                 (1.0_dp/Real((num_symm*(N_in+1)),dp))*factor*((real(g(1)*conjg(g(1)),dp) + &
                                 real(g(2)*conjg(g(2)),dp) + real(g(3)*conjg(g(3)),dp))/3.0_dp)
                         end do
                      end do
                   end if
                elseif(index(optics_geom,'tensor')>0)then
                   if (num_symm==0) then 
                      do N2=1,3  
                         g(N2) = band_gradient(n_eigen,n_eigen2,N2,N,N_spin)
                      end do
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,1) = factor*real(g(1)*conjg(g(1)),dp)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,2) = factor*real(g(2)*conjg(g(2)),dp)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,3) = factor*real(g(3)*conjg(g(3)),dp)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,4) = factor*real(g(1)*conjg(g(2)),dp)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,5) = factor*real(g(1)*conjg(g(3)),dp)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,6) = factor*real(g(2)*conjg(g(3)),dp)
                   else
                      do N2=1,num_symm
                         do N3=1,1+N_in
                            g = 0.0_dp
                            do i=1,3  
                               do j=1,3
                                  g(i) = g(i) + ((-1.0_dp)**(N3+1))*(crystal_symmetry_operations(i,j,N2)&
                                       *band_gradient(n_eigen,n_eigen2,j,N,N_spin))
                               end do
                            end do
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,1) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,1) + factor*real(g(1)*conjg(g(1)),dp)/Real((num_symm*(N_in+1)),dp)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,2) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,2) + factor*real(g(2)*conjg(g(2)),dp)/Real((num_symm*(N_in+1)),dp)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,3) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,3) + factor*real(g(3)*conjg(g(3)),dp)/Real((num_symm*(N_in+1)),dp)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,4) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,4) + factor*real(g(1)*conjg(g(2)),dp)/Real((num_symm*(N_in+1)),dp)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,5) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,5) + factor*real(g(1)*conjg(g(3)),dp)/Real((num_symm*(N_in+1)),dp)
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,6) =  matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,6) + factor*real(g(2)*conjg(g(3)),dp)/Real((num_symm*(N_in+1)),dp)

                         end do
                      end do
                   end if
                endif
             end do
          end do
       end do
    end do
  end subroutine make_weights

  !***************************************************************
  subroutine calc_epsilon_2
    !***************************************************************
    ! This subroutine calculates epsilon_2

    use od_constants, only : dp, pi
    use od_cell, only : nkpoints, cell_volume 
    use od_electronic, only : nspins, electrons_per_state, nbands
    use od_jdos_utils, only : E,jdos_nbins
    use od_parameters, only : optics_intraband, optics_drude_broadening

    integer :: N_energy
    integer :: N
    integer :: N_spin
    integer :: N2

    real(kind=dp) ::dE
    real(kind=dp) :: x
    real(kind=dp) :: epsilon2_const
 
    dE = E(2)-E(1)
    epsilon2_const = (e_charge*pi*1E-20)/(cell_volume*1E-30*epsilon_0)

    if(optics_intraband)then
       allocate(intra(N_geom))
       do N=1,N_geom
          do N_spin=1,nspins
             intra(N_geom) = intra(N_geom) + weighted_dos_at_e(N_spin,N_geom)
          enddo
       enddo
       intra = intra*e_charge/(cell_volume*1E-10*epsilon_0)
    endif

    if(.not. optics_intraband) then 
       allocate(epsilon(jdos_nbins,2,N_geom,1))
    else
       allocate(epsilon(jdos_nbins,2,N_geom,3))
    endif
    epsilon=0.0_dp

    do N2=1,N_geom
       do N_spin=1,nspins                        ! Loop over spins
          do N_energy=2,jdos_nbins
             epsilon(N_energy,2,N2,1) = epsilon(N_energy,2,N2,1) + &
                  epsilon2_const*weighted_jdos(N_energy,N_spin,N2)
             if (optics_intraband) then 
                epsilon(N_energy,2,N2,2) = epsilon(N_energy,2,N2,2) + &
                     ((intra(N2)*(e_charge**2)*hbar*optics_drude_broadening)/(((E(N_energy)*e_charge)**2)+((optics_drude_broadening*hbar)**2))) 
                epsilon(N_energy,2,N2,3) = epsilon(N_energy,2,N2,3) + &
                     epsilon(N_energy,2,N2,2) + epsilon(N_energy,2,N2,1)*E(N_energy)*e_charge
             end if
          end do
       end do
    end do

    ! Sum rule 
    if (N_geom==1) then 
       x = 0.0_dp
       do N=2,jdos_nbins   !! don't include 0eV as it makes in intraband case blow up (and should be zero otherwise)
          if(.not. optics_intraband)then
             x = x+((N*(dE**2)*epsilon(N,2,1,1))/(hbar**2))
          else
             x = x+((N*(dE**2)*epsilon(N,2,1,3))/((hbar**2)*E(N)*e_charge))
          endif
       end do
       N_eff = (x*e_mass*cell_volume*1E-30*epsilon_0*2)/(pi) 
    end if

  end subroutine calc_epsilon_2

  !***************************************************************
  subroutine calc_epsilon_1
    !***************************************************************
    ! This subroutine uses kramers kronig to calculate epsilon_1 

    use od_constants, only : dp, pi
    use od_jdos_utils, only : E, jdos_nbins
    use od_parameters, only : optics_intraband, optics_drude_broadening

    integer :: N_energy
    integer :: N_energy2
    integer :: N2
    real(kind=dp),allocatable, dimension(:) :: q
    real(kind=dp) :: energy1
    real(kind=dp) :: energy2
    real(kind=dp) :: dE

    dE=E(2)-E(1)
    if(.not. optics_intraband) then
       allocate(q(1))
    else 
       allocate(q(3))
    endif

    do N2=1,N_geom
       do N_energy=1,jdos_nbins
          q=0.0_dp
          do N_energy2=1,jdos_nbins                 
             if (N_energy2.ne.N_energy) then
                energy1 = E(N_energy)  
                energy2 = E(N_energy2)
                q(1)=q(1)+(((energy2*epsilon(N_energy2,2,N2,1))/((energy2**2)-(energy1**2)))*dE)
                if(optics_intraband)then 
                   q(2)=q(2)+((dE*epsilon(N_energy2,2,N2,2))/(((energy2**2)-(energy1**2))*e_charge))
                   q(3)=q(3)+((dE*epsilon(N_energy2,2,N2,3))/(((energy2**2)-(energy1**2))*e_charge))
                endif
             end if
          end do
          epsilon(N_energy,1,N2,1)=((2.0_dp/pi)*q(1))+1.0_dp
          if(optics_intraband) then 
!             epsilon(N_energy,1,N2,2)=((2.0_dp/pi)*q(2))+1.0_dp  !! old KK method 
             epsilon(N_energy,1,N2,2)=1.0_dp-(intra(N_geom)/((E(N_energy)**2)+(((optics_drude_broadening*hbar)/e_charge)**2)))
!             epsilon(N_energy,1,N2,3)=((2.0_dp/pi)*q(3))+1.0_dp  !! old KK method
             epsilon(N_energy,1,N2,3)=epsilon(N_energy,1,N2,1)+epsilon(N_energy,1,N2,2)-1.0_dp
         endif
       end do
    end do

  end subroutine calc_epsilon_1

  !***************************************************************
  subroutine calc_loss_fn
    !***************************************************************
    ! This subroutine calculates the loss function and the sum rules

    use od_constants, only : dp, cmplx_i, pi
    use od_jdos_utils, only : E, jdos_nbins
    use od_cell, only : cell_volume
    use od_parameters, only : optics_intraband

    complex(kind=dp) :: g 
    integer :: N_energy
    real(kind=dp) :: x
    real(kind=dp) :: dE

    if(.not. optics_intraband) allocate(loss_fn(jdos_nbins,1)) 
    if(optics_intraband) allocate(loss_fn(jdos_nbins,3)) 
    loss_fn=0.0_dp

    dE=E(2)-E(1)
    g = (0.0_dp,0.0_dp)

    do N_energy=1,jdos_nbins
       g = epsilon(N_energy,1,1,1)+(cmplx_i*epsilon(N_energy,2,1,1))
       loss_fn(N_energy,1)=-1*aimag(1.0_dp/g)
       if(optics_intraband) then 
          g = epsilon(N_energy,1,1,2)+(cmplx_i*epsilon(N_energy,2,1,2)/(E(N_energy)*e_charge))
          loss_fn(N_energy,2)=-1*aimag(1.0_dp/g)  
          g = epsilon(N_energy,1,1,3)+(cmplx_i*epsilon(N_energy,2,1,3)/(E(N_energy)*e_charge))
          loss_fn(N_energy,3)=-1*aimag(1.0_dp/g)
       endif
    end do
    loss_fn(1,2)=0.0_dp  ! gets rid of the NaN from dividing by zero in the loop above
    loss_fn(1,3)=0.0_dp

    ! Sum rule 1
    x = 0.0_dp
    do N_energy=2,jdos_nbins
       if(.not. optics_intraband)then 
          x = x+(N_energy*(dE**2)*loss_fn(N_energy,1))
       else
          x = x+(N_energy*(dE**2)*loss_fn(N_energy,3))
       endif
    end do
    N_eff2 = x*(e_mass*cell_volume*1E-30*epsilon_0*2)/(pi*(hbar**2))

    ! Sum rule 2
    x = 0
    do N_energy=2,jdos_nbins
       if(.not. optics_intraband)then 
          x = x+(loss_fn(N_energy,1)/N_energy)
       else
          x = x+(loss_fn(N_energy,3)/N_energy)
       endif
    end do
    N_eff3 = x

  end subroutine calc_loss_fn

  !***************************************************************
  subroutine calc_conduct
    !***************************************************************
    ! This subroutine calculates the conductivity

    use od_jdos_utils, only : jdos_nbins, E
    use od_parameters, only : optics_intraband 

    integer :: N_energy

    allocate(conduct(1:jdos_nbins,2))  
    conduct=0.0_dp

    if(.not. optics_intraband)then
       do N_energy=1,jdos_nbins
          conduct(N_energy,1)=(E(N_energy)*e_charge/hbar)*epsilon_0*epsilon(N_energy,2,1,1)
       end do
       do N_energy=1,jdos_nbins
          conduct(N_energy,2)=(E(N_energy)*e_charge/hbar)*epsilon_0*(1.0_dp-epsilon(N_energy,1,1,1))
       end do
    else
       do N_energy=1,jdos_nbins
          conduct(N_energy,1)=epsilon_0*epsilon(N_energy,2,1,3)/hbar
       end do
       do N_energy=1,jdos_nbins
          conduct(N_energy,2)=(E(N_energy)*e_charge/hbar)*epsilon_0*(1.0_dp-epsilon(N_energy,1,1,3))
       end do
    end if

  end subroutine calc_conduct

  !***************************************************************
  subroutine calc_refract
    !***************************************************************
    ! This subroutine calculates the refractive index

    use od_jdos_utils, only : jdos_nbins, E
    use od_parameters, only : optics_intraband

    integer :: N_energy

    allocate(refract(jdos_nbins,2)) 
    refract=0.0_dp

    if(.not. optics_intraband)then
       do N_energy=1,jdos_nbins
          refract(N_energy,1)=(0.5_dp*((((epsilon(N_energy,1,1,1)**2)+&
               &(epsilon(N_energy,2,1,1)**2))**0.5_dp)+epsilon(N_energy,1,1,1)))**(0.5_dp)
       end do
       do N_energy=1,jdos_nbins
          refract(N_energy,2)=(0.5_dp*((((epsilon(N_energy,1,1,1)**2)+&
               &(epsilon(N_energy,2,1,1)**2))**0.5_dp)-epsilon(N_energy,1,1,1)))**(0.5_dp)
       end do
    else 
       do N_energy=1,jdos_nbins
          refract(N_energy,1)=(0.5_dp*((((epsilon(N_energy,1,1,3)**2)+&
               &((epsilon(N_energy,2,1,3)/(E(N_energy)*e_charge))**2))**0.5_dp)+epsilon(N_energy,1,1,1)))**(0.5_dp)
       end do
       do N_energy=1,jdos_nbins
          refract(N_energy,2)=(0.5_dp*((((epsilon(N_energy,1,1,1)**2)+&
               &((epsilon(N_energy,2,1,3)/(E(N_energy)*e_charge))**2))**0.5_dp)-epsilon(N_energy,1,1,1)))**(0.5_dp)
       end do

    end if

  end subroutine calc_refract

  !***************************************************************
  subroutine calc_absorp
    !***************************************************************
    ! This subroutine calculates the absorption coefficient

    use od_jdos_utils, only : jdos_nbins, E

    integer :: N_energy

    allocate(absorp(jdos_nbins)) 
    absorp=0.0_dp

    do N_energy=1,jdos_nbins
       absorp(N_energy)=2*refract(N_energy,2)*E(N_energy)*e_charge/(hbar*c_speed)
    end do

  end subroutine calc_absorp

  !***************************************************************
  subroutine calc_reflect
    !***************************************************************
    ! This subroutine calculates the reflection coefficient

    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy

    allocate(reflect(jdos_nbins)) 
    reflect=0.0_dp

    do N_energy=1,jdos_nbins
       reflect(N_energy)=(((refract(N_energy,1)-1)**2)+(refract(N_energy,2)**2))/&
            &(((refract(N_energy,1)+1)**2)+(refract(N_energy,2)**2))
    end do

  end subroutine calc_reflect

  !***************************************************************
  subroutine write_epsilon
    !***************************************************************
    ! This subroutine writes out the dielectric function

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, &
         &output_format, fixed, adaptive, linear, optics_intraband
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only: E, jdos_nbins
    use od_io, only : seedname, io_file_unit, stdout

    integer :: N,N2
    real(kind=dp) ::dE
    integer :: epsilon_unit

    type(graph_labels) :: label

    label%name="epsilon"
    label%title="Dielectric Function"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="Real"
    label%legend_b="Imaginary"

    dE=E(2)-E(1)

    ! Open the output file
    epsilon_unit = io_file_unit()
    open(unit=epsilon_unit,action='write',file=trim(seedname)//'.epsilon')

    ! Write into the output file
    write(epsilon_unit,*)'#*********************************************'
    write(epsilon_unit,*)'#            Dielectric function                 '
    write(epsilon_unit,*)'#*********************************************'
    write(epsilon_unit,*)'#'
    write(epsilon_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(epsilon_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(epsilon_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(epsilon_unit,*)'# Number of bands:',nbands
    write(epsilon_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(epsilon_unit,*)'#'
    write(epsilon_unit,'(1x,a,f10.6,1x,a,f10.6,1x,a)')'# Dielectric function calculated to',jdos_max_energy,'eV in',dE,'eV steps'
    write(epsilon_unit,*)'#'
    write(epsilon_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(epsilon_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(epsilon_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(epsilon_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(epsilon_unit,*)'#'
    if(optics_intraband)then 
       write(epsilon_unit,*)'# Calculation includes intraband term'
       if(fixed) write(epsilon_unit,*)'# DOS at Ef:', dos_at_e(1,:)
       if(adaptive) write(epsilon_unit,*)'# DOS at Ef:', dos_at_e(2,:)
       if(linear) write(epsilon_unit,*)'# DOS at Ef:', dos_at_e(3,:)
       write(epsilon_unit,*)'# Plasmon energy:', (intra(1)**0.5)
    endif
    if (N_geom==1) then
       write(epsilon_unit,*)'# Result of sum rule: Neff(E) =  ',N_eff
       write(epsilon_unit,*)'#'
       if(.not. optics_intraband)then  
          do N=1,jdos_nbins
             write(epsilon_unit,*)E(N),epsilon(N,1,1,1),epsilon(N,2,1,1)
          enddo
       else 
          do N2=1,3
             write(epsilon_unit,*)''
             write(epsilon_unit,*)''         
             do N=1,jdos_nbins
                write(epsilon_unit,*)E(N),epsilon(N,1,1,N2),epsilon(N,2,1,N2)
             enddo
          enddo
       endif

       if(trim(output_format)=="xmgrace") then
          call write_optics_xmgrace(label,E,epsilon(:,1,1,1),epsilon(:,2,1,1))
       elseif(trim(output_format)=="gnuplot") then 
          write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
       else
          write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
       endif

    end if
    if (index(optics_geom,'tensor')>0) then
       do N2=1,N_geom  
          write(epsilon_unit,*)''
          write(epsilon_unit,*)''
          do N=1,jdos_nbins
             write(epsilon_unit,*)E(N),epsilon(N,1,N2,1),epsilon(N,2,N2,1)
          end do

          label%name="epsilon"//trim(achar(N2))
          if(trim(output_format)=="xmgrace") then
             call write_optics_xmgrace(label,E,epsilon(:,1,N2,1),epsilon(:,2,N2,1))
          elseif(trim(output_format)=="gnuplot") then 
             write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
          else
             write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
          endif
       end do
    end if

    ! Close output file 
    close(unit=epsilon_unit)  

  end subroutine write_epsilon

  !***************************************************************
  subroutine write_loss_fn
    !***************************************************************
    ! This subroutine writes out the loss function

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, output_format, &
         optics_intraband
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only: seedname, io_file_unit,stdout

    integer :: N, N2 
    integer :: loss_fn_unit

    type(graph_labels) :: label


    label%name="loss_fn"
    label%title="Loss Function"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="loss_fn"

    ! Open the output file
    loss_fn_unit = io_file_unit()
    open(unit=loss_fn_unit,action='write',file=trim(seedname)//'.loss_fn')

    ! Write into the output file
    write(loss_fn_unit,*)'#*********************************************'
    write(loss_fn_unit,*)'#               Loss function                 '
    write(loss_fn_unit,*)'#*********************************************'
    write(loss_fn_unit,*)'#'
    write(loss_fn_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(loss_fn_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(loss_fn_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(loss_fn_unit,*)'# No of bands:',nbands
    write(loss_fn_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(loss_fn_unit,*)'#'
    write(loss_fn_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(loss_fn_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(loss_fn_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(loss_fn_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(loss_fn_unit,*)'#'
    write(loss_fn_unit,*)'# Result of first sum rule: Neff(E) = ',N_eff2
    write(loss_fn_unit,*)'# Result of second sum rule (pi/2 = 1.570796327):',N_eff3
    write(loss_fn_unit,*)'#' 
    if(.not. optics_intraband) then 
       do N=1,jdos_nbins
          write(loss_fn_unit,*)E(N),loss_fn(N,1)
       end do
    else
       do N2=1,3
          do N=1,jdos_nbins 
             write(loss_fn_unit,*)E(N),loss_fn(N,N2)
          enddo
       enddo
    endif


    ! Close output file 
    close(unit=loss_fn_unit)  

    if(trim(output_format)=="xmgrace") then
       !    call write_optics_xmgrace(label,E,loss_fn) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif(trim(output_format)=="gnuplot") then 
       write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
    else
       write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
    endif

  end subroutine write_loss_fn

  !***************************************************************
  subroutine write_conduct
    !***************************************************************
    ! This subroutine writes out the conductivity

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, output_format
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit,stdout

    integer :: N 
    integer :: conduct_unit

    type(graph_labels) :: label


    label%name="conductivity"
    label%title="Conductivity"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="Real"
    label%legend_b="Imaginary"


    ! Open the output file
    conduct_unit = io_file_unit()
    open(unit=conduct_unit,action='write',file=trim(seedname)//'.conductivity')

    ! Write into the output file
    write(conduct_unit,*)'#*********************************************'
    write(conduct_unit,*)'#               Conductivity                '
    write(conduct_unit,*)'#*********************************************'
    write(conduct_unit,*)'#'
    write(conduct_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(conduct_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(conduct_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(conduct_unit,*)'# No of bands:',nbands
    write(conduct_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(conduct_unit,*)'#'
    write(conduct_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(conduct_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(conduct_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(conduct_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(conduct_unit,*)'#'
    do N=1,jdos_nbins
       write(conduct_unit,*)E(N),conduct(N,1),conduct(N,2)
    end do

    ! Close output file
    close(unit=conduct_unit)  

    if(trim(output_format)=="xmgrace") then
       call write_optics_xmgrace(label,E,conduct(:,1),conduct(:,2))
    elseif(trim(output_format)=="gnuplot") then 
       write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
    else
       write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
    endif

  end subroutine write_conduct

  !***************************************************************
  subroutine write_refract
    !***************************************************************
    ! This subroutine writes out the refractive index

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, output_format
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit, stdout

    integer :: N 
    integer :: refract_unit

    type(graph_labels) :: label


    label%name="refractive_index"
    label%title="Refractive Index"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="Real"
    label%legend_b="Imaginary"

    ! Open the output file
    refract_unit = io_file_unit()
    open(unit=refract_unit,action='write',file=trim(seedname)//'.refractive_index')

    ! Write into the output file
    write(refract_unit,*)'#*********************************************'
    write(refract_unit,*)'#             Refractive index                 '
    write(refract_unit,*)'#*********************************************'
    write(refract_unit,*)'#'
    write(refract_unit,*)'# N=n+ik'
    write(refract_unit,*)'#'
    write(refract_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(refract_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(refract_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(refract_unit,*)'# No of bands:',nbands
    write(refract_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(refract_unit,*)'#'
    write(refract_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(refract_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(refract_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(refract_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(refract_unit,*)'#'    
    do N=1,jdos_nbins
       write(refract_unit,*)E(N),refract(N,1),refract(N,2)
    end do

    ! Close output file
    close(unit=refract_unit)  

    if(trim(output_format)=="xmgrace") then
       call write_optics_xmgrace(label,E,refract(:,1),refract(:,2))
    elseif(trim(output_format)=="gnuplot") then 
       write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
    else
       write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
    endif

  end subroutine write_refract

  !***************************************************************
  subroutine write_absorp
    !***************************************************************
    ! This subroutine writes out the absorption coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, output_format
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit, stdout

    integer :: N 
    integer :: absorp_unit

    type(graph_labels) :: label


    label%name="absorption"
    label%title="Absorption Coefficient"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="A"

    ! Open the output file
    absorp_unit = io_file_unit()
    open(unit=absorp_unit,action='write',file=trim(seedname)//'.absorption')

    ! Write into the output file
    write(absorp_unit,*)'#*********************************************'
    write(absorp_unit,*)'#             Absorption coefficent                 '
    write(absorp_unit,*)'#*********************************************'
    write(absorp_unit,*)'#'
    write(absorp_unit,*)'#'
    write(absorp_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(absorp_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(absorp_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(absorp_unit,*)'# No of bands:',nbands
    write(absorp_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(absorp_unit,*)'#'
    write(absorp_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(absorp_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(absorp_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(absorp_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(absorp_unit,*)'#'    
    do N=1,jdos_nbins
       write(absorp_unit,*)E(N),absorp(N)
    end do

    ! Close output file
    close(unit=absorp_unit)  

    if(trim(output_format)=="xmgrace") then
       call write_optics_xmgrace(label,E,absorp)
    elseif(trim(output_format)=="gnuplot") then 
       write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
    else
       write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
    endif

  end subroutine write_absorp

  !***************************************************************
  subroutine write_reflect
    !***************************************************************
    ! This subroutine writes out the reflection coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op, output_format
    use od_electronic, only: nbands, num_electrons, nspins
    use od_io, only : seedname, io_file_unit, stdout
    use od_jdos_utils, only : jdos_nbins, E

    integer :: N
    integer :: reflect_unit
    type(graph_labels) :: label


    label%name="reflection"
    label%title="Reflection Coefficient"
    label%x_label="Energy (eV)"
    label%y_label=""
    label%legend_a="R"

    ! Open the output file
    reflect_unit = io_file_unit()
    open(unit=reflect_unit,action='write',file=trim(seedname)//'.reflection')

    ! Write into the output file
    write(reflect_unit,*)'#*********************************************'
    write(reflect_unit,*)'#           Reflection coefficient                '
    write(reflect_unit,*)'#*********************************************'
    write(reflect_unit,*)'#'
    write(reflect_unit,*)'# N=n+ik'
    write(reflect_unit,*)'#'
    write(reflect_unit,*)'# Number of k-points: ',nkpoints 
    if(nspins==1)then
       write(reflect_unit,*)'# Number of electrons:',num_electrons(1)
    else
       write(reflect_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    end if
    write(reflect_unit,*)'# No of bands:',nbands
    write(reflect_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(reflect_unit,*)'#'
    write(reflect_unit,*)'# optics_geom:  ',optics_geom
    if (index(optics_geom,'polar')>0) then 
       write(reflect_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
       write(reflect_unit,*)'# q_weight:',q_weight
    end if
    if (scissor_op>0) then 
       write(reflect_unit,'(1x,a,f10.3,f10.3,f10.3)')'# Scissor operator:', scissor_op
    end if
    write(reflect_unit,*)'#'    
    do N=1,jdos_nbins
       write(reflect_unit,*)E(N),reflect(N)
    end do

    ! Close output file 
    close(unit=reflect_unit)  

    if(trim(output_format)=="xmgrace") then
       call write_optics_xmgrace(label,E,reflect)
    elseif(trim(output_format)=="gnuplot") then 
       write(stdout,*)  " WARNING: GNUPLOT output not yet available, continuing..."
    else
       write(stdout,*)  " WARNING: Unknown output format requested, continuing..."
    endif


  end subroutine write_reflect

  !=============================================================================== 
  subroutine write_optics_xmgrace(label,E,column1,column2)
    !=============================================================================== 
    use xmgrace_utils
    use od_io,         only : io_file_unit,io_error,seedname 
    implicit none 

    type(graph_labels),intent(in) :: label

    !  type :: graph_labels
    !     character(10) :: name
    !     character(20) :: title
    !     character(20) :: x_label
    !     character(20) :: y_label
    !     character(20) :: legend_a
    !     character(20) :: legend_b
    !  end type graph_labels

    real(dp),  intent(in) :: E(:)
    real(dp),  intent(in)  :: column1(:)
    real(dp),  optional, intent(in) :: column2(:)

    real(dp) :: min_x, max_x, min_y, max_y

    integer :: batch_file,ierr, array_lengths

    batch_file=io_file_unit()
    open(unit=batch_file,file=trim(seedname)//'.'//trim(label%name)//'.agr',iostat=ierr)
    if(ierr.ne.0) call io_error(" ERROR: Cannot open xmgrace batch file in dos: write_dos_xmgrace")

    min_x=minval(E)
    max_x=maxval(E)


    if(present(column2)) then
       min_y=min(minval(column1), minval(column2))
    else
       max_y=minval(column1)
    endif


    if(present(column2)) then
       max_y=max(maxval(column1), maxval(column2))
    else
       max_y=maxval(column1)
    endif


    call  xmgu_setup(batch_file)
    call  xmgu_legend(batch_file)
    call  xmgu_title(batch_file, min_x, max_x, min_y, max_y, trim(label%title))
    call  xmgu_subtitle(batch_file,"Generated by OptaDOS")

    call  xmgu_axis(batch_file,"x",trim(label%x_label))
    call  xmgu_axis(batch_file,"y",trim(label%y_label))

    if(present(column2)) then
       call xmgu_data_header(batch_file,0,1,trim(label%legend_a))
       call xmgu_data_header(batch_file,1,2,trim(label%legend_b))
       call xmgu_data(batch_file,0,E(:),column1(:))
       call xmgu_data(batch_file,1,E(:),column2(:))
    else
       call xmgu_data_header(batch_file,0,1,trim(label%legend_a))
       call xmgu_data(batch_file,0,E(:),column1)  
    endif

  end subroutine write_optics_xmgrace

endmodule od_optics
