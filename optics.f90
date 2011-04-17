!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
module od_optics

  implicit none
  private
  public :: optics_calculate

  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp) :: q_weight

  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_jdos

  real(kind=dp),allocatable, dimension(:,:) :: epsilon
  real(kind=dp),allocatable, dimension(:,:) :: conduct
  real(kind=dp),allocatable, dimension(:,:) :: refract
  real(kind=dp),allocatable, dimension(:,:) :: loss_fn
  real(kind=dp),allocatable, dimension(:,:) :: absorp
  real(kind=dp),allocatable, dimension(:,:) :: reflect
    
  real(kind=dp) :: N_eff
  real(kind=dp) :: N_eff2
  real(kind=dp) :: N_eff3

  real(kind=dp), parameter :: epsilon_0 = 8.8541878176E-12  ! need to put in correct number
  real(kind=dp), parameter :: e_charge =  1.60217646E-19 ! need to put in correct number 
  real(kind=dp), parameter :: e_mass =  9.10938E-31 ! need to put in correct number
  real(kind=dp), parameter :: hbar =  1.054571628E-34 ! need to put in correct number 
  real(kind=dp), parameter :: c_speed =  299792458 ! need to put in correct number

contains

  subroutine optics_calculate
    !
    !  Program to calculate optics properties
    !

    use od_constants, only : dp
    use od_electronic, only : band_gradient, elec_read_band_gradient
    use od_cell, only : cell_volume
    use od_jdos_utils, only : jdos_utils_calculate
    use od_comms, only : on_root

    !! Get information from .bands file 
    !! Will need to get information on the kpoints    

    ! Get information from .cst_ome file 
    call elec_read_band_gradient

    !! Get symmetry will also go here 

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
       call calc_conduct
       call calc_refract
       call calc_loss_fn
       call calc_absorp
       call calc_reflect
       
       ! Write everything out
       call write_epsilon
       call write_conduct
       call write_refract
       call write_loss_fn
       call write_absorp
       call write_reflect
       
    endif

  end subroutine optics_calculate

  ! Subroutines go here 

  !***************************************************************
  subroutine make_weights
    !***************************************************************
    use od_constants, only : dp
    use od_electronic, only : nbands, nspins, band_gradient, num_electrons, electrons_per_state
    use od_cell, only : nkpoints, cell_volume, num_kpoints_on_node
    use od_parameters, only : optics_geom, optics_qdir
    use od_io, only : io_error
    use od_comms, only : my_node_id

    real(kind=dp), dimension(3) :: qdir 
    integer :: N
    integer :: N_spin
    integer :: n_eigen
    integer :: n_eigen2
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp) :: g

    num_occ = 0.0_dp
    DO N_spin=1,nspins
       num_occ(N_spin) = num_electrons(N_spin)
    ENDDO

    if (electrons_per_state==2) then 
       num_occ(1) = num_occ(1)/2.0_dp
    endif

    allocate(matrix_weights(nbands,nbands,num_kpoints_on_node(my_node_id),nspins))
    matrix_weights=0.0_dp

    qdir=optics_qdir
    q_weight=((qdir(1)**2.0_dp)+(qdir(2)**2.0_dp)+(qdir(3)**2.0_dp))**0.5_dp
    if(q_weight<0.001_dp)&
         call io_error("Error:  please check optics_qdir, norm close to zero")

    DO N=1,num_kpoints_on_node(my_node_id)                      ! Loop over kpoints
       DO N_spin=1,nspins                                    ! Loop over spins
          DO n_eigen=1,NINT(num_occ(N_spin))                 ! Loop over occupied states 
             DO n_eigen2=(NINT(num_occ(N_spin))+1),nbands    ! Loop over unoccupied states
                g = (((qdir(1)*band_gradient(n_eigen,n_eigen2,1,N,N_spin))+ &
                     (qdir(2)*band_gradient(n_eigen,n_eigen2,2,N,N_spin))+ &
                     (qdir(3)*band_gradient(n_eigen,n_eigen2,3,N,N_spin)))/q_weight)
                matrix_weights(n_eigen,n_eigen2,N,N_spin) = real(g*conjg(g),dp)
                !           matrix_weights(n_eigen,n_eigen2,N,N_spin) = 1.0_dp  ! for JDOS 
             END DO
          END DO
       END DO
    END DO
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

    real(kind=dp) ::dE
    real(kind=dp) :: x
    real(kind=dp) :: epsilon2_const

    dE = E(2)-E(1)
    epsilon2_const = (e_charge*pi*1E-20)/(cell_volume*1E-30*epsilon_0)

    allocate(epsilon(jdos_nbins,3))
    epsilon=0.0_dp

    epsilon(1,3) = 0.0_dp ! set epsilon_2=0 at 0eV

!    DO N=1,nkpoints                               ! Loop over kpoints
       DO N_spin=1,nspins                        ! Loop over spins
          DO N_energy=2,jdos_nbins
             epsilon(N_energy,3)=epsilon(N_energy,3) + epsilon2_const*(1.0_dp/(E(N_energy)**2))*weighted_jdos(N_energy,N_spin,1)
             !           epsilon(N_energy,3)=epsilon(N_energy,3) + weighted_jdos(N_energy,N_spin,N) !output the jDOS 
          END DO
       END DO
!    END DO

    ! Sum rule 
    x = 0.0_dp
    DO N=1,jdos_nbins
       x = x+((N*(dE**2)*epsilon(N,3))/(hbar**2))
    END DO
    N_eff = (x*e_mass*cell_volume*1E-30*epsilon_0*2)/(pi) 

  end subroutine calc_epsilon_2

  !***************************************************************
  subroutine calc_epsilon_1
    !***************************************************************
    ! This subroutine uses kramers kronig to calculate epsilon_1 

    use od_constants, only : dp, pi
    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy
    integer :: N_energy2
    real(kind=dp) :: q 
    real(kind=dp) :: energy1
    real(kind=dp) :: energy2
    real(kind=dp) :: dE

    dE=E(2)-E(1)

    DO N_energy=1,jdos_nbins
       q=0.0_dp
       DO N_energy2=1,jdos_nbins                  ! NOTE This starts at 2 at the moment
          IF (N_energy2.ne.N_energy) then
             energy1 = E(N_energy)  
             energy2 = E(N_energy2)
             q=q+(((energy2*epsilon(N_energy2,3))/((energy2**2)-(energy1**2)))*dE)
          END IF
       END DO
       epsilon(N_energy,2)=((2.0_dp/pi)*q)+1.0_dp
    END DO

  end subroutine calc_epsilon_1

  !***************************************************************
  subroutine calc_loss_fn
    !***************************************************************
    ! This subroutine calculates the loss function and the sum rules

    use od_constants, only : dp, cmplx_i, pi
    use od_jdos_utils, only : E,jdos_nbins
    use od_cell, only : cell_volume

    complex(kind=dp) :: g 
    integer :: N_energy
    real(kind=dp) :: x
    real(kind=dp) :: dE

    allocate(loss_fn(1:jdos_nbins,2)) 
    loss_fn=0.0_dp

    dE=E(2)-E(1)

    DO N_energy=1,jdos_nbins
       loss_fn(N_energy,1)=E(N_energy)
    END DO

    g = (0.0_dp,0.0_dp)

    DO N_energy=1,jdos_nbins
       g = epsilon(N_energy,2)+(cmplx_i*epsilon(N_energy,3))
       loss_fn(N_energy,2)=-1*AIMAG(1.0_dp/g)
    END DO

    ! Sum rule 1
    x = 0.0_dp
    DO N_energy=1,jdos_nbins
       x = x+(N_energy*(dE**2)*loss_fn(N_energy,2))
    END DO
    N_eff2 = x*(e_mass*cell_volume*1E-30*epsilon_0*2)/(pi*(hbar**2))

    ! Sum rule 2
    x = 0
    DO N_energy=1,jdos_nbins
       x = x+(loss_fn(N_energy,2)/N_energy)
    END DO 
    N_eff3 = x

  end subroutine calc_loss_fn

  !***************************************************************
  subroutine calc_conduct
    !***************************************************************
    ! This subroutine calculates the conductivity

    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy

    allocate(conduct(1:jdos_nbins,3)) 
    conduct=0.0_dp

    DO N_energy=1,jdos_nbins
       conduct(N_energy,1)=E(N_energy)
    END DO

    DO N_energy=1,jdos_nbins
       conduct(N_energy,2)=(conduct(N_energy,1)*e_charge/hbar)*epsilon_0*epsilon(N_energy,3)
    END DO

    DO N_energy=1,jdos_nbins
       conduct(N_energy,3)=(conduct(N_energy,1)*e_charge/hbar)*epsilon_0*(1.0_dp-epsilon(N_energy,2))
    END DO

  end subroutine calc_conduct

  !***************************************************************
  subroutine calc_refract
    !***************************************************************
    ! This subroutine calculates the refractive index

    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy

    allocate(refract(jdos_nbins,3)) 
    refract=0.0_dp

    DO N_energy=1,jdos_nbins
       refract(N_energy,1)=E(N_energy)
    END DO

    DO N_energy=1,jdos_nbins
       refract(N_energy,2)=(0.5_dp*((((epsilon(N_energy,2)**2.0_dp)+&
            &(epsilon(N_energy,3)**2.0_dp))**0.5_dp)+epsilon(N_energy,2)))**(0.5_dp)
    END DO

    DO N_energy=1,jdos_nbins
       refract(N_energy,3)=(0.5_dp*((((epsilon(N_energy,2)**2.0_dp)+&
            &(epsilon(N_energy,3)**2.0_dp))**0.5_dp)-epsilon(N_energy,2)))**(0.5_dp)
    END DO

  end subroutine calc_refract

  !***************************************************************
  subroutine calc_absorp
    !***************************************************************
    ! This subroutine calculates the absorption coefficient

    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy

    allocate(absorp(jdos_nbins,2)) 
    absorp=0.0_dp

    DO N_energy=1,jdos_nbins
       absorp(N_energy,1)=E(N_energy)
    END DO

    DO N_energy=1,jdos_nbins
       absorp(N_energy,2)=2*refract(N_energy,3)*absorp(N_energy,1)*e_charge/(hbar*c_speed)
    END DO

  end subroutine calc_absorp

  !***************************************************************
  subroutine calc_reflect
    !***************************************************************
    ! This subroutine calculates the reflection coefficient

    use od_jdos_utils, only : E,jdos_nbins

    integer :: N_energy

    allocate(reflect(jdos_nbins,2)) 
    reflect=0.0_dp

    DO N_energy=1,jdos_nbins
       reflect(N_energy,1)=E(N_energy)
    END DO

    DO N_energy=1,jdos_nbins
       reflect(N_energy,2)=(((refract(N_energy,2)-1)**2)+(refract(N_energy,3)**2))/&
            &(((refract(N_energy,2)+1)**2)+(refract(N_energy,3)**2))
    END DO

  end subroutine calc_reflect

  !***************************************************************
  subroutine write_epsilon
    !***************************************************************
    ! This subroutine writes out the dielectric function

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only: E, jdos_nbins
    use od_io, only : seedname, io_file_unit

    integer :: N
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
    IF(nspins==1)THEN
       write(epsilon_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(epsilon_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(epsilon_unit,*)'# Number of bands:',nbands
    write(epsilon_unit,*)'#Volume of the unit cell (Ang^3):',cell_volume
    write(epsilon_unit,*)'#'
    write(epsilon_unit,'(1x,a,f10.6,1x,a,f10.6,1x,a)')'# Dielectric function calculated to',jdos_max_energy,'eV in',dE,'eV steps'
    write(epsilon_unit,*)'#'
    write(epsilon_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(epsilon_unit,*)'# q_weight:',q_weight
    write(epsilon_unit,*)'#'
    write(epsilon_unit,*)'# Result of sum rule: Neff(E) =  ',N_eff
    write(epsilon_unit,*)'#'    
    DO N=1,jdos_nbins
       write(epsilon_unit,*)E(N),epsilon(N,2),epsilon(N,3)
    END DO

    ! Close output file 
    close(unit=epsilon_unit)  

  end subroutine write_epsilon

  !***************************************************************
  subroutine write_loss_fn
    !***************************************************************
    ! This subroutine writes out the loss function

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins
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
    IF(nspins==1)THEN
       write(loss_fn_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(loss_fn_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(loss_fn_unit,*)'# No of bands:',nbands
    write(loss_fn_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(loss_fn_unit,*)'#'
    write(loss_fn_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(loss_fn_unit,*)'# q_weight',q_weight
    write(loss_fn_unit,*)'#'
    write(loss_fn_unit,*)'# Result of first sum rule: Neff(E) = ',N_eff2
    write(loss_fn_unit,*)'# Result of second sum rule (pi/2 = 1.570796327):',N_eff3
    write(loss_fn_unit,*)'#'    
    DO N=1,jdos_nbins
       write(loss_fn_unit,*)loss_fn(N,1),loss_fn(N,2)
    END DO

    ! Close output file 
    close(unit=loss_fn_unit)  

  end subroutine write_loss_fn

  !***************************************************************
  subroutine write_conduct
    !***************************************************************
    ! This subroutine writes out the conductivity

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins
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
    IF(nspins==1)THEN
       write(conduct_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(conduct_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(conduct_unit,*)'# No of bands:',nbands
    write(conduct_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(conduct_unit,*)'#'
    write(conduct_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(conduct_unit,*)'# q_weight',q_weight
    write(conduct_unit,*)'#'
    DO N=1,jdos_nbins
       write(conduct_unit,*)conduct(N,1),conduct(N,2),conduct(N,3)
    END DO

    ! Close output file
    close(unit=conduct_unit)  

  end subroutine write_conduct

  !***************************************************************
  subroutine write_refract
    !***************************************************************
    ! This subroutine writes out the refractive index

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins
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
    IF(nspins==1)THEN
       write(refract_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(refract_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(refract_unit,*)'# No of bands:',nbands
    write(refract_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(refract_unit,*)'#'
    write(refract_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(refract_unit,*)'# q_weight',q_weight
    write(refract_unit,*)'#'    
    DO N=1,jdos_nbins
       write(refract_unit,*)refract(N,1),refract(N,2),refract(N,3)
    END DO

    ! Close output file
    close(unit=refract_unit)  

  end subroutine write_refract

  !***************************************************************
  subroutine write_absorp
    !***************************************************************
    ! This subroutine writes out the absorption coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_jdos_utils, only : jdos_nbins
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
    IF(nspins==1)THEN
       write(absorp_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(absorp_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(absorp_unit,*)'# No of bands:',nbands
    write(absorp_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(absorp_unit,*)'#'
    write(absorp_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(absorp_unit,*)'# q_weight',q_weight
    write(absorp_unit,*)'#'    
    DO N=1,jdos_nbins
       write(absorp_unit,*)absorp(N,1),absorp(N,2)
    END DO

    ! Close output file
    close(unit=absorp_unit)  

  end subroutine write_absorp

  !***************************************************************
  subroutine write_reflect
    !***************************************************************
    ! This subroutine writes out the reflection coefficient

    use od_cell, only : nkpoints, cell_volume
    use od_parameters, only : optics_geom, optics_qdir,jdos_max_energy
    use od_electronic, only: nbands, num_electrons, nspins
    use od_io, only : seedname, io_file_unit
    use od_jdos_utils, only : jdos_nbins

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
    IF(nspins==1)THEN
       write(reflect_unit,*)'# Number of electrons:',num_electrons(1)
    ELSE
       write(reflect_unit,*)'# Number of electrons:',num_electrons(1),num_electrons(2)
    END IF
    write(reflect_unit,*)'# No of bands:',nbands
    write(reflect_unit,*)'# Volume of the unit cell (Ang^3):',cell_volume
    write(reflect_unit,*)'#'
    write(reflect_unit,'(1x,a,f10.3,f10.3,f10.3)')'# q-vector', optics_qdir(1),optics_qdir(2),optics_qdir(3)
    write(reflect_unit,*)'# q_weight',q_weight
    write(reflect_unit,*)'#'    
    DO N=1,jdos_nbins
       write(reflect_unit,*)reflect(N,1),reflect(N,2)
    END DO

    ! Close output file 
    close(unit=reflect_unit)  

  end subroutine write_reflect

endmodule od_optics
