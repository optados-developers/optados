!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
module od_optics

  implicit none
  private
  public :: optics_calculate

  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp) :: q_weight 

  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_jdos

  real(kind=dp),allocatable, dimension(:,:,:) :: epsilon
  real(kind=dp),allocatable, dimension(:,:) :: conduct
  real(kind=dp),allocatable, dimension(:,:) :: refract
  real(kind=dp),allocatable, dimension(:) :: loss_fn
  real(kind=dp),allocatable, dimension(:) :: absorp
  real(kind=dp),allocatable, dimension(:) :: reflect

  real(kind=dp) :: N_eff
  real(kind=dp) :: N_eff2
  real(kind=dp) :: N_eff3
  integer :: N_geom

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
    use od_electronic, only : band_gradient, elec_read_band_gradient
    use od_cell, only : cell_volume
    use od_jdos_utils, only : jdos_utils_calculate
    use od_comms, only : on_root
    use od_parameters, only : optics_geom

    ! Get information from .cst_ome file 
    call elec_read_band_gradient

    ! Form matrix element
    call make_weights

    ! Send matrix element to jDOS routine and get weighted jDOS back
    call jdos_utils_calculate(matrix_weights, weighted_jdos)

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
          do n_eigen=1,nint(num_occ(N_spin))                 ! Loop over occupied states 
             do n_eigen2=(nint(num_occ(N_spin))+1),nbands    ! Loop over unoccupied states
                if(index(optics_geom,'unpolar')>0)then
                   if (num_symm==0) then 
                      g(1) = (((qdir1(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir1(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir1(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight1)
                      g(2) = (((qdir2(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir2(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir2(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = &
                           0.5_dp*(real(g(1)*conjg(g(1)),dp)+real(g(2)*conjg(g(2)),dp))/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)&
                           +scissor_op)**2)
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
                                 (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)/&
                                 ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                                 scissor_op)**2) 
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
                                 (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)/&
                                 ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                                 scissor_op)**2) 
                         end do
                      end do
                   end if
                elseif(index(optics_geom,'polar')>0)then
                   if (num_symm==0) then 
                      g(1) = (((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                           (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                           (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = real(g(1)*conjg(g(1)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
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
                                 (1.0_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)/&
                                 ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                                 scissor_op)**2) 
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
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,N_geom) = (real(g(1)*conjg(g(1)),dp)+&
                           real(g(2)*conjg(g(2)),dp) + real(g(3)*conjg(g(3)),dp))/&
                           (((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)*3.0_dp)
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
                                 (1.0_dp/Real((num_symm*(N_in+1)),dp))*(real(g(1)*conjg(g(1)),dp) + &
                                 real(g(2)*conjg(g(2)),dp) + real(g(3)*conjg(g(3)),dp))/&
                                 (((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                                 scissor_op)**2)*3.0_dp)
                         end do
                      end do
                   end if
                elseif(index(optics_geom,'tensor')>0)then
                   if (num_symm==0) then 
                      do N2=1,3  
                         g(N2) = band_gradient(n_eigen,n_eigen2,N2,N,N_spin)
                      end do
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,1) = real(g(1)*conjg(g(1)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,2) = real(g(2)*conjg(g(2)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,3) = real(g(3)*conjg(g(3)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,4) = real(g(1)*conjg(g(2)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,5) = real(g(1)*conjg(g(3)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
                      matrix_weights(n_eigen,n_eigen2,N,N_spin,6) = real(g(2)*conjg(g(3)),dp)/&
                           ((band_energy(n_eigen2,N_spin,N)-band_energy(n_eigen,N_spin,N)+&
                           scissor_op)**2)
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
                                 N,N_spin,1) + real(g(1)*conjg(g(1)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,2) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,2) + real(g(2)*conjg(g(2)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,3) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,3) + real(g(3)*conjg(g(3)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,4) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,4) + real(g(1)*conjg(g(2)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,5) = matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,5) + real(g(1)*conjg(g(3)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))
                            matrix_weights(n_eigen,n_eigen2,N,N_spin,6) =  matrix_weights(n_eigen,n_eigen2,&
                                 N,N_spin,6) + real(g(2)*conjg(g(3)),dp)/(((band_energy(n_eigen2,N_spin,N)&
                                 -band_energy(n_eigen,N_spin,N)+scissor_op)**2)*Real((num_symm*(N_in+1)),dp))     
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
    use od_electronic, only : nspins, electrons_per_state
    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy
    integer :: N
    integer :: N_spin
    integer :: N2

    real(kind=dp) ::dE
    real(kind=dp) :: x
    real(kind=dp) :: epsilon2_const

    dE = E(2)-E(1)
    epsilon2_const = (e_charge*pi*1E-20)/(cell_volume*1E-30*epsilon_0)

    allocate(epsilon(jdos_nbins,2,N_geom))
    epsilon=0.0_dp

    do N2=1,N_geom
       epsilon(1,2,N2) = 0.0_dp ! set epsilon_2=0 at 0eV
    end do

    do N2=1,N_geom
       do N_spin=1,nspins                        ! Loop over spins
          do N_energy=2,jdos_nbins
             epsilon(N_energy,2,N2) = epsilon(N_energy,2,N2) + &
                  epsilon2_const*weighted_jdos(N_energy,N_spin,N2)
          end do
       end do
    end do

    ! Sum rule 
    if (N_geom==1) then  
       x = 0.0_dp
       do N=1,jdos_nbins
          x = x+((N*(dE**2)*epsilon(N,2,1))/(hbar**2))
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

    integer :: N_energy
    integer :: N_energy2
    integer :: N2
    real(kind=dp) :: q 
    real(kind=dp) :: energy1
    real(kind=dp) :: energy2
    real(kind=dp) :: dE

    dE=E(2)-E(1)

    do N2=1,N_geom
       do N_energy=1,jdos_nbins
          q=0.0_dp
          do N_energy2=1,jdos_nbins                 
             if (N_energy2.ne.N_energy) then
                energy1 = E(N_energy)  
                energy2 = E(N_energy2)
                q=q+(((energy2*epsilon(N_energy2,2,N2))/((energy2**2)-(energy1**2)))*dE)
             end if
          end do
          epsilon(N_energy,1,N2)=((2.0_dp/pi)*q)+1.0_dp
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

    complex(kind=dp) :: g 
    integer :: N_energy
    real(kind=dp) :: x
    real(kind=dp) :: dE

    allocate(loss_fn(jdos_nbins)) 
    loss_fn=0.0_dp

    dE=E(2)-E(1)
    g = (0.0_dp,0.0_dp)

    do N_energy=1,jdos_nbins
       g = epsilon(N_energy,1,1)+(cmplx_i*epsilon(N_energy,2,1))
       loss_fn(N_energy)=-1*aimag(1.0_dp/g)
    end do

    ! Sum rule 1
    x = 0.0_dp
    do N_energy=1,jdos_nbins
       x = x+(N_energy*(dE**2)*loss_fn(N_energy))
    end do
    N_eff2 = x*(e_mass*cell_volume*1E-30*epsilon_0*2)/(pi*(hbar**2))

    ! Sum rule 2
    x = 0
    do N_energy=1,jdos_nbins
       x = x+(loss_fn(N_energy)/N_energy)
    end do
    N_eff3 = x

  end subroutine calc_loss_fn

  !***************************************************************
  subroutine calc_conduct
    !***************************************************************
    ! This subroutine calculates the conductivity

    use od_jdos_utils, only : jdos_nbins, E

    integer :: N_energy

    allocate(conduct(1:jdos_nbins,2))  
    conduct=0.0_dp

    do N_energy=1,jdos_nbins
       conduct(N_energy,1)=(E(N_energy)*e_charge/hbar)*epsilon_0*epsilon(N_energy,2,1)
    end do

    do N_energy=1,jdos_nbins
       conduct(N_energy,2)=(E(N_energy)*e_charge/hbar)*epsilon_0*(1.0_dp-epsilon(N_energy,1,1))
    end do

  end subroutine calc_conduct

  !***************************************************************
  subroutine calc_refract
    !***************************************************************
    ! This subroutine calculates the refractive index

    use od_jdos_utils, only : jdos_nbins

    integer :: N_energy

    allocate(refract(jdos_nbins,2)) 
    refract=0.0_dp

    do N_energy=1,jdos_nbins
       refract(N_energy,1)=(0.5_dp*((((epsilon(N_energy,1,1)**2)+&
            &(epsilon(N_energy,2,1)**2))**0.5_dp)+epsilon(N_energy,1,1)))**(0.5_dp)
    end do

    do N_energy=1,jdos_nbins
       refract(N_energy,2)=(0.5_dp*((((epsilon(N_energy,1,1)**2)+&
            &(epsilon(N_energy,2,1)**2))**0.5_dp)-epsilon(N_energy,1,1)))**(0.5_dp)
    end do

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
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only: E, jdos_nbins
    use od_io, only : seedname, io_file_unit

    integer :: N,N2
    real(kind=dp) ::dE
    integer :: epsilon_unit

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
    if (N_geom==1) then
       write(epsilon_unit,*)'# Result of sum rule: Neff(E) =  ',N_eff
       write(epsilon_unit,*)'#'   
       do N=1,jdos_nbins
          write(epsilon_unit,*)E(N),epsilon(N,1,1),epsilon(N,2,1)
       end do
    end if
    if (index(optics_geom,'tensor')>0) then
       do N2=1,N_geom  
          write(epsilon_unit,*)''
          write(epsilon_unit,*)''
          do N=1,jdos_nbins
             write(epsilon_unit,*)E(N),epsilon(N,1,N2),epsilon(N,2,N2)
          end do
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
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only: seedname, io_file_unit

    integer :: N 
    integer :: loss_fn_unit

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
    do N=1,jdos_nbins
       write(loss_fn_unit,*)E(N),loss_fn(N)
    end do

    ! Close output file 
    close(unit=loss_fn_unit)  

  end subroutine write_loss_fn

  !***************************************************************
  subroutine write_conduct
    !***************************************************************
    ! This subroutine writes out the conductivity

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit

    integer :: N 
    integer :: conduct_unit

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

  end subroutine write_conduct

  !***************************************************************
  subroutine write_refract
    !***************************************************************
    ! This subroutine writes out the refractive index

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit

    integer :: N 
    integer :: refract_unit

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

  end subroutine write_refract

  !***************************************************************
  subroutine write_absorp
    !***************************************************************
    ! This subroutine writes out the absorption coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : seedname, io_file_unit

    integer :: N 
    integer :: absorp_unit

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

  end subroutine write_absorp

  !***************************************************************
  subroutine write_reflect
    !***************************************************************
    ! This subroutine writes out the reflection coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy, scissor_op
    use od_electronic, only: nbands, num_electrons, nspins
    use od_io, only : seedname, io_file_unit
    use od_jdos_utils, only : jdos_nbins, E

    integer :: N
    integer :: reflect_unit

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

  end subroutine write_reflect

endmodule od_optics
