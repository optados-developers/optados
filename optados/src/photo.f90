!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
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
!=============================================================================== 
module od_photo

  use od_constants, only : dp

  implicit none
  private
  public :: photo_calculate


  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: pdos_weights_atoms
  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: pdos_weights_atoms_tmp
  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: projected_matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: optical_matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:,:,:) :: foptical_matrix_weights
  real(kind=dp),allocatable, public, dimension(:,:,:) :: weighted_jdos
  real(kind=dp),allocatable, public, dimension(:,:) :: absorp_layer
  real(kind=dp),allocatable, dimension(:,:,:) :: pdos_weights_k_band
  real(kind=dp),allocatable, dimension(:,:,:) :: imfp_val
  real(kind=dp),allocatable, dimension(:,:,:,:) :: electron_esc
  real(kind=dp), dimension(:,:), allocatable :: I_layer
  real(kind=dp), dimension(:,:), allocatable :: absorption_layer
  real(kind=dp), dimension(:), allocatable :: total_absorption
  real(kind=dp), dimension(:), allocatable :: total_transmittance
  real(kind=dp), allocatable, public, save :: E(:)
  real(kind=dp), save                   :: delta_bins
  integer,save :: jdos_nbins
  real(kind=dp),allocatable, dimension(:) :: intra
  real(kind=dp),allocatable, public, dimension(:,:) :: weighted_dos_at_e_photo
  real(kind=dp),allocatable, dimension(:,:,:,:) :: epsilon
  real(kind=dp),allocatable, dimension(:,:) :: epsilon_sum 
  real(kind=dp),allocatable, dimension(:,:,:) :: refract_photo
  real(kind=dp),allocatable, dimension(:,:)  :: reflect_photo
  real(kind=dp),allocatable, dimension(:,:) :: absorp_photo

  real(kind=dp),allocatable, dimension(:,:) :: refract
  real(kind=dp),allocatable, dimension(:)  :: reflect
  real(kind=dp),allocatable, dimension(:) :: absorp

  real(kind=dp), dimension(:), allocatable :: thickness_atom
  real(kind=dp), dimension(:,:), allocatable :: new_atoms_coordinates
  real(kind=dp),allocatable, dimension(:,:,:) :: phi_arpes
  real(kind=dp),allocatable, dimension(:,:,:) :: theta_arpes
  real(kind=dp),allocatable, dimension(:,:,:) :: E_kinetic
  real(kind=dp),allocatable, dimension(:,:,:) :: E_transverse
  real(kind=dp),allocatable, dimension(:,:,:) :: bulk_prob
  real(kind=dp),allocatable, dimension(:) :: t_energy
  real(kind=dp),allocatable, dimension(:,:) :: weighted_qe
  real(kind=dp), allocatable, dimension(:,:,:,:,:) :: weighted_temp
  integer :: max_energy
  real(kind=dp),allocatable, dimension(:,:,:,:) :: qe_osm
  real(kind=dp),allocatable, dimension(:,:,:,:,:) :: qe_tsm
  integer, dimension(:), allocatable :: atom_order
  integer, dimension(:), allocatable :: atoms_per_layer
  real (kind=dp), dimension(3) :: t_cart
  real (kind=dp) :: work_function_eff
  real (kind=dp) :: evacuum
  real (kind=dp) :: evacuum_eff
  real (kind=dp) :: total_field_emission
  real(kind=dp),allocatable, dimension(:,:,:) :: field_emission
  integer,allocatable, dimension(:) :: layer
  integer :: N_geom
  integer :: max_atoms
  integer :: max_doubling_atom
  integer :: max_layer
  real(kind=dp) :: e_fermi
  real(kind=dp) :: q_weight 
  contains

  subroutine photo_calculate
    !
    !  Program to calculate the photoemission
    !

   use od_electronic, only : elec_dealloc_optical, elec_pdos_read, &
        optical_mat, elec_read_optical_mat, nbands, nspins, &
        efermi, efermi_set, elec_dealloc_optical,elec_read_foptical_mat
   use od_cell, only : cell_volume,num_kpoints_on_node
   use od_jdos_utils, only : jdos_utils_calculate,setup_energy_scale
   use od_comms, only : on_root, my_node_id
   use od_parameters, only : optics_geom, adaptive, linear, fixed, optics_intraband, &
        optics_drude_broadening, work_function,photo_model,&
        elec_field
    use od_dos_utils, only :  dos_utils_set_efermi, dos_utils_calculate_at_e
    use od_io, only : stdout,io_error 
    use od_pdos, only : pdos_calculate




    if(on_root) then
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+                             Photoemission Calculation                      +'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '|                                                                            |'
    endif

    if(.not.efermi_set) call dos_utils_set_efermi


    !Identify layers

    call calc_layers

    !THIS PART COMES FROM THE PDOS MODULE
    ! read in the pdos weights
    call elec_pdos_read
   
    call make_pdos_weights_atoms
     
    call elec_read_optical_mat

    ! Calculate the optical properties of the slab
    call calc_photo_optics

    call calc_absorp_layer

    !Calculate the photoemission angles theta/phi and transverse energy
    call calc_angle

    !Calculate the electron escape lenght
    call calc_electron_esc

    call bulk_emission

    !Electric field and field emission
    if(elec_field.gt.0.0_dp)then
        call effect_wf
    else
        evacuum_eff=efermi+work_function
        work_function_eff = work_function
    end if

    !Calculate the QE
    if (index(photo_model,'3step')>0) then !Three-step-model
       call calc_three_step_model
    end if
    if (index(photo_model,'1step')>0) then !One-step-model
       call elec_read_foptical_mat !Read the one-step matrix elements
       call make_foptical_weights !Calculate the one-step optical matrix
       call calc_one_step_model !Calculate QE
    end if

    call weighted_mean_te !Weight the contribution of each electron
                         !to the transverse energy spread according to their QE
    !Broaden ouputs using a gaussian function
    call binding_energy_spread
    !Write either a binding energy output with after Gaussian broadening
    call write_qe_output_files

    !Deallocate everything
    call photo_deallocate    

    write(stdout,*)'End of execution'


   end subroutine photo_calculate


    !***************************************************************
    subroutine calc_layers            
    !***************************************************************
    !This subroutine identifies the layer of each atom
        use od_constants, only : dp
        use od_cell, only :  num_atoms, atoms_pos_cart_photo, atoms_label_tmp             
        use od_io, only : stdout, io_error

     integer :: atom_1,atom_2,i,index,temp,first,ierr,atom

    allocate(atom_order(num_atoms),stat=ierr)
    if(ierr/=0) call io_error('Error: calc_layers - allocation of atom_order failed')
    
    allocate(layer(num_atoms),stat=ierr)
    if(ierr/=0) call io_error('Error: calc_layers - allocation of layer failed')

    do i=1,num_atoms
        atom_order(i)=i
    end do
     
    !SORTING ALGORITHM
    do atom_1=1,num_atoms-1
        first = atom_order(atom_1)
        do atom_2=atom_1+1,num_atoms
           index = atom_1
              if(atoms_pos_cart_photo(3,atom_order(atom_2)).gt.atoms_pos_cart_photo(3,first)) then
                  first = atom_order(atom_2)
                  index = atom_2     
              end if

              if(index /= atom_1) then
                  temp = atom_order(atom_1)
                  atom_order(atom_1) = atom_order(index)
                  atom_order(index) = temp
              end if
        end do
    end do


    !DEFINE THE LAYER FOR EACH ATOM
    i=1
    layer(1) = 1
    do atom=2,num_atoms
        if(( trim(atoms_label_tmp(atom_order(atom))).ne.trim(atoms_label_tmp(atom_order(atom-1)))).or.& 
        (ABS(atoms_pos_cart_photo(3,atom_order(atom))-atoms_pos_cart_photo(3,atom_order(atom-1))).gt.0.50))then
            i = i+1
        end if
        layer(atom)=i
    end do

    write(stdout,'(1x,a78)') '+------------------------------- Atomic Order  ------------------------------+'
    write(stdout,'(1x,a78)') '| Atom |  Atom Order  |   Layer   |          Atomic Position (Angs)          |'

    do atom=1,num_atoms
       write(stdout,*) "|  ",trim(atoms_label_tmp(atom_order(atom))),atom_order(atom),&
                   layer(atom),'              ',&
                   atoms_pos_cart_photo(3,atom_order(atom)),"      |"
    end do

    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'

    !CALCULATE THE MAX LAYER (HALF SLAB)
    max_layer=((layer(num_atoms)+1)/2)


    !CALCULATE THE MAX ATOM (HALF SLAB)
    max_atoms=0
    do atom=1,num_atoms
        if(layer(atom).le.((layer(num_atoms)+1)/2))then
            max_atoms = max_atoms+1
        end if
    end do



    write(stdout,*)'|  Max number of atoms:', max_atoms,'   Max  number of layers:', max_layer,'   |'

    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
   
    allocate(atoms_per_layer(max_layer),stat=ierr)
    if(ierr/=0) call io_error('Error: calc_layers - allocation of atoms_per_layer failed')

    !CALCULATE HOW MANY ATOMS PER LAYER
    atoms_per_layer=1
    do atom=2,max_atoms
        if(layer(atom).eq.layer(atom-1))then
            atoms_per_layer(layer(atom))=atoms_per_layer(layer(atom))+1
        end if
    end do

  end subroutine calc_layers            

  !***************************************************************
  subroutine make_pdos_weights_atoms
    !***************************************************************
    !This subroutine is equivalent to pdos_merge of pdos.F90, but only for atoms
    use od_electronic, only :  pdos_orbital, pdos_weights, pdos_mwab, nspins, &
    num_electrons, efermi, band_energy
    use od_cell, only : num_kpoints_on_node, num_atoms
    use od_comms, only : my_node_id
    use od_io, only : io_error

    integer :: N,N_spin,n_eigen,np,ierr,atom,i,j

    allocate(pdos_weights_atoms(num_atoms,pdos_mwab%nbands,num_kpoints_on_node(my_node_id),nspins),stat=ierr)
        if(ierr/=0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')

    allocate(pdos_weights_atoms_tmp(num_atoms,pdos_mwab%nbands,num_kpoints_on_node(my_node_id),nspins),stat=ierr)
        if(ierr/=0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_atoms failed')


    allocate(pdos_weights_k_band(pdos_mwab%nbands,num_kpoints_on_node(my_node_id),nspins),stat=ierr)
        if(ierr/=0) call io_error('Error: make_pdos_weights_atoms - allocation of pdos_weights_k_band failed')


        pdos_weights_atoms=0.0_dp
        pdos_weights_k_band=0.0_dp

    do N=1,num_kpoints_on_node(my_node_id) 
        do N_spin=1,nspins                 
            do n_eigen=1,pdos_mwab%nbands  
                i=1
                do np=1,pdos_mwab%norbitals
                    if(np.gt.1)then
                        if(pdos_orbital%rank_in_species(np).ne.pdos_orbital%rank_in_species(np-1))then
                            i=i+1
                        end if
                    end if
                    pdos_weights_atoms(i,n_eigen,N,N_spin)=&
                    pdos_weights_atoms(i,n_eigen,N,N_spin)+&
                    pdos_weights(np,n_eigen,N,N_spin)
                end do
            end do
        end do
    end do
    
    do N=1,num_kpoints_on_node(my_node_id) 
        do N_spin=1,nspins                 
            do n_eigen=1,pdos_mwab%nbands  
               do atom=1,max_atoms
                  if(pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin).lt.0.0_dp)then
                    pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin)=0.0_dp
                  end if
                    pdos_weights_k_band(n_eigen,N,N_spin)=pdos_weights_k_band(n_eigen,N,N_spin)+&
                    pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin)
                end do
            end do
        end do
    end do

  end subroutine make_pdos_weights_atoms

    !***************************************************************
    subroutine calc_photo_optics
    !***************************************************************

    use od_optics, only : make_weights,calc_epsilon_2,calc_epsilon_1,&
        calc_refract,calc_absorp,calc_reflect,epsilon,refract,absorp,reflect,&
        intra
    use od_io, only : stdout, io_error
    use od_electronic, only : optical_mat, elec_read_optical_mat, nbands,nspins,&
          efermi, efermi_set, elec_dealloc_optical,&
          band_gradient,elec_read_band_gradient,&
          nbands, nspins, band_energy
    use od_cell, only : cell_volume, num_kpoints_on_node, kpoint_r,&
          num_atoms, num_kpoints_on_node
    use od_jdos_utils, only : jdos_utils_calculate,jdos_nbins, E,&
          setup_energy_scale
    use od_comms, only : on_root, my_node_id
    use od_parameters, only : optics_geom, adaptive, linear, fixed, optics_intraband, &
          optics_drude_broadening,slab_volume
    use od_dos_utils, only : dos_utils_calculate_at_e
    use od_constants, only : epsilon_0,e_charge

    real(kind=dp),allocatable, dimension(:,:,:,:) :: dos_matrix_weights
    real(kind=dp),allocatable, dimension(:,:) :: weighted_dos_at_e
    real(kind=dp),allocatable, dimension(:,:) :: weighted_dos_at_e_photo
    real(kind=dp),allocatable, dimension(:,:) :: dos_at_e
    integer :: N, N2,N_spin,n_eigen,n_eigen2,atom,ierr,N_energy



    allocate(absorp_photo(jdos_nbins,max_atoms))
    allocate(reflect_photo(jdos_nbins,max_atoms))



    call make_weights(matrix_weights)

    do atom=1,max_atoms                           ! Loop over atoms 

    allocate(projected_matrix_weights(nbands,nbands,num_kpoints_on_node(my_node_id),nspins,N_geom))
    if(ierr/=0) call io_error('Error: make_photo_weights - allocation of matrix_weights failed')

    projected_matrix_weights = 0.0_dp
      do N2 = 1, N_geom
        do N=1,num_kpoints_on_node(my_node_id)    ! Loop over kpoints
            do N_spin=1,nspins                    ! Loop over spins
                do n_eigen=1,nbands               ! Loop over state 1 
                    do n_eigen2=n_eigen,nbands    ! Loop over state 2 
                    if(band_energy(n_eigen,N_spin,N)>efermi.and.n_eigen/=n_eigen2)cycle
                    if(band_energy(n_eigen2,N_spin,N)<efermi.and.n_eigen/=n_eigen2)cycle
                        projected_matrix_weights(n_eigen,n_eigen2,N,N_spin,N2) = &
                        matrix_weights(n_eigen,n_eigen2,N,N_spin,N2)*&
                        (pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin)/pdos_weights_k_band(n_eigen,N,N_spin))
                    end do                        ! Loop over state 2 
                end do                            ! Loop over state 1
            end do                                ! Loop over spins
        end do                                    ! Loop over kpoints
     end do

      ! Send matrix element to jDOS routine and get weighted jDOS back
      call jdos_utils_calculate(projected_matrix_weights, weighted_jdos)

      if(allocated(projected_matrix_weights)) then
         deallocate(projected_matrix_weights,stat=ierr)
         if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate projected_matrix_weights')
      end if

    if (optics_intraband) then
      allocate (dos_matrix_weights(size(matrix_weights, 5), nbands,num_kpoints_on_node(my_node_id), nspins))
      allocate (dos_at_e(3, nspins))
      allocate (weighted_dos_at_e(nspins, size(matrix_weights, 5)))
      weighted_dos_at_e = 0.0_dp
      do N = 1, size(matrix_weights, 5)
        do N2 = 1, nbands
          dos_matrix_weights(N, N2, :, :) = matrix_weights(N2, N2, :, :, N)
        enddo
      enddo
      call dos_utils_calculate_at_e(efermi, dos_at_e, dos_matrix_weights, weighted_dos_at_e)
    endif

    if (on_root) then
      ! Calculate epsilon_2
      call calc_epsilon_2(weighted_jdos,weighted_dos_at_e)

      ! Calculate epsilon_1
      call calc_epsilon_1

      ! Calculate other optical properties
        call calc_refract
        call calc_absorp
        call calc_reflect

     endif

       absorp_photo(N_energy,atom) = absorp(N_energy)

       reflect_photo(N_energy,atom) = reflect(N_energy)

         deallocate(dos_matrix_weights)
         deallocate(dos_at_e)
         deallocate(weighted_dos_at_e)
         deallocate(weighted_jdos)
         deallocate(epsilon)
         deallocate(refract)
         deallocate(absorp)
         deallocate(reflect)
         deallocate(E)
     if(allocated(intra)) then
         deallocate(intra)
     end if
    end do                                        ! Loop over atoms

  end subroutine calc_photo_optics

  !===============================================================================
  subroutine make_foptical_weights    
  !===============================================================================
  ! This subroutine calclualtes te optical matrix elements for the one step
  ! photoemission model.
  ! Victor Chang, 7th February 2020
  !===============================================================================

    use od_constants, only : dp
    use od_electronic, only : nbands, nspins, optical_mat, num_electrons, &
         electrons_per_state, band_energy, efermi,foptical_mat
    use od_cell, only : nkpoints, cell_volume, num_kpoints_on_node, cell_get_symmetry, &
         num_crystal_symmetry_operations, crystal_symmetry_operations, num_atoms
    use od_parameters, only : optics_geom, optics_qdir, legacy_file_format, scissor_op, devel_flag,photon_energy
    use od_io, only : io_error
    use od_comms, only : my_node_id

    real(kind=dp), dimension(3) :: qdir 
    real(kind=dp), dimension(3) :: qdir1 
    real(kind=dp), dimension(3) :: qdir2
    real(kind=dp) :: q_weight1 
    real(kind=dp) :: q_weight2 
    integer :: N, i, j, N_2
    integer :: N_in
    integer :: N_spin, N_spin_2
    integer :: N2, N3
    integer :: n_eigen, n_eigen_2
    integer :: n_eigen2, n_eigen2_2
    integer :: num_symm
    integer :: ierr
    real(kind=dp), dimension(2) :: num_occ
    complex(kind=dp), dimension(3) :: g
    real(kind=dp) :: factor
    real(kind=dp) :: test



    if (.not.legacy_file_format.and.index(devel_flag,'old_filename')>0) then 
        num_symm=0
        call cell_get_symmetry
    end if 
    num_symm = num_crystal_symmetry_operations

    num_occ = 0.0_dp
    do N_spin=1,nspins
        num_occ(N_spin) = num_electrons(N_spin)
    enddo

    if (electrons_per_state==2) then 
        num_occ(1) = num_occ(1)/2.0_dp
    endif

       N_geom=1
 
    allocate(foptical_matrix_weights(nbands+1,nbands+1,num_kpoints_on_node(my_node_id),nspins,N_geom))
    if(ierr/=0) call io_error('Error: make_optical_weights - allocation of optical_matrix_weights failed')
    foptical_matrix_weights=0.0_dp

    if (index(optics_geom,'polar')>0) then 
        qdir=optics_qdir
        q_weight=((qdir(1)**2)+(qdir(2)**2)+(qdir(3)**2))**0.5_dp
        if(q_weight<0.001_dp)&
        call io_error("Error:  please check optics_qdir, norm close to zero")
    endif

    if (index(optics_geom,'unpolar')>0) then 
        !TO CHANGE WHEN THE light_direction IS CORRECTED
        !optics_qdir(:)=t_cart(:)
        if(optics_qdir(3).lt.1E-06)then 
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
                 factor = 1.0_dp/(photon_energy**2) 
                        if(index(optics_geom,'unpolar')>0)then
                            if (num_symm==0) then 
                                g(1) = (((qdir1(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                (qdir1(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                (qdir1(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight1)
                                g(2) = (((qdir2(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                (qdir2(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                (qdir2(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight2)
                                foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) = &
                                0.5_dp*factor*(real(g(1)*conjg(g(1)),dp)+real(g(2)*conjg(g(2)),dp))
                            else ! begin unpolar symmetric
                                do N2=1,num_symm
                                    do N3=1,1+N_in
                                        do i=1,3
                                            qdir(i)=0.0_dp
                                            do j=1,3
                                                qdir(i) = qdir(i) + ((-1.0_dp)**(N3+1))*&
                                                (crystal_symmetry_operations(j,i,N2)*qdir1(j))
                                            end do
                                        end do
                                        g(1)=(((qdir(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                        (qdir(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                        (qdir(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight1)
                                        foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) =  &
                                        foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) + &
                                        (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)*factor
                                        g(1) = 0.0_dp
                                            do i=1,3 ! if I include an extra variable I can merge this and the last do loops
                                                qdir(i)=0.0_dp
                                                do j=1,3
                                                    qdir(i) = qdir(i) + ((-1.0_dp)**(N3+1))*&
                                                    (crystal_symmetry_operations(j,i,N2)*qdir2(j))
                                                end do
                                            end do
                                        g(1)=(((qdir(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                        (qdir(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                        (qdir(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight2)
                                        foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) =  &
                                        foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) + &
                                        (0.5_dp/Real((num_symm*(N_in+1)),dp))*real(g(1)*conjg(g(1)),dp)*factor
                                    end do
                                end do
                            end if !end unpolar symmetric
                        elseif(index(optics_geom,'polar')>0)then
                            if (num_symm==0) then 
                                g(1) = (((qdir(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                (qdir(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                (qdir(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight)
                                foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) = factor*real(g(1)*conjg(g(1)),dp)
                            else !begin polar symmetric
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
                                    g(1)=(((qdir(1)*foptical_mat(n_eigen,nbands+1,1,N,N_spin))+ &
                                    (qdir(2)*foptical_mat(n_eigen,nbands+1,2,N,N_spin))+ &
                                    (qdir(3)*foptical_mat(n_eigen,nbands+1,3,N,N_spin)))/q_weight)
                                    foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) =  &
                                    foptical_matrix_weights(n_eigen,nbands+1,N,N_spin,N_geom) + &
                                    (1.0_dp/Real((num_symm*(N_in+1)),dp))*factor*real(g(1)*conjg(g(1)),dp)
                                    end do
                                end do
                            end if !end polar symmetric
                    end if ! end photo_geom
            end do       ! Loop over state 1 
        end do           ! Loop over spins
    end do               ! Loop over kpoints

  end subroutine make_foptical_weights

  !***************************************************************
  subroutine calc_absorp_layer
    !***************************************************************
    ! This subroutine calculates the absorption coefficient

    use od_cell, only : num_atoms, atoms_pos_cart_photo, num_species
    use od_jdos_utils, only : jdos_nbins, E
    use od_parameters, only : photon_energy,&
    jdos_spacing,surface_area
    use od_io, only : stdout, io_error

    real(kind=dp), dimension(:), allocatable :: light_path      
    real(kind=dp), dimension(:,:), allocatable :: attenuation_layer
    real(kind=dp) :: transmittance
    integer :: atom,i,N_energy,ierr,first_atom_second_l,last_atom_secondlast_l
    real(kind=dp) :: I_0

    allocate(thickness_atom(max_atoms))
    if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of thickness_atom failed')

    thickness_atom=0.0_dp

        allocate(light_path(max_atoms))
        if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of light_path failed')


        light_path=0.0_dp

        allocate(I_layer(jdos_nbins,max_layer))
        if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of I_layer failed')


    allocate(attenuation_layer(jdos_nbins,max_atoms))
        if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of attenuation_layer failed')

        allocate(absorption_layer(jdos_nbins,max_atoms))
        if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of absorption_layer  failed')


    allocate(total_absorption(jdos_nbins))
    if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of total_absorption failed')
    allocate(total_transmittance(jdos_nbins))
    if(ierr/=0) call io_error('Error: calc_absorp_layer - allocation of total_transmittance failed')


    if(max_layer.lt.2)then
        thickness_atom=1.5
    else

        do atom=2,max_atoms
            if(layer(atom).gt.1)then
                thickness_atom(1)=((atoms_pos_cart_photo(3,atom_order(1))-atoms_pos_cart_photo(3,atom_order(atom)))/2)*2
                first_atom_second_l=atom
                exit
            end if
        end do
    
        do i=2,first_atom_second_l-1
            thickness_atom(i)=((atoms_pos_cart_photo(3,atom_order(i))-&
            atoms_pos_cart_photo(3,atom_order(first_atom_second_l)))/2)*2
        end do
            
        do i=1,max_atoms
            if(layer(max_atoms-i).lt.layer(max_atoms))then
                thickness_atom(max_atoms)=(ABS(atoms_pos_cart_photo(3,atom_order(max_atoms))-&
                atoms_pos_cart_photo(3,atom_order(max_atoms-i)))/2)*2
                last_atom_secondlast_l=max_atoms-i
                exit
            end if
        end do
    
        do i=last_atom_secondlast_l+1,max_atoms-1
            thickness_atom(i)=(ABS(atoms_pos_cart_photo(3,atom_order(i))-&
            atoms_pos_cart_photo(3,atom_order(last_atom_secondlast_l)))/2)*2
        end do
    
    
        do atom=first_atom_second_l,last_atom_secondlast_l
            thickness_atom(atom)=abs((atoms_pos_cart_photo(3,atom_order(sum(atoms_per_layer(1:layer(atom)-1))))&
            -atoms_pos_cart_photo(3,atom_order(atom)))/2)+&
            abs((atoms_pos_cart_photo(3,atom_order(atom))-&
            atoms_pos_cart_photo(3,atom_order(sum(atoms_per_layer(1:layer(atom)))+1)))/2)  
      end do
    end if

        do atom=1,max_atoms
            light_path(atom)=thickness_atom(atom)
        end do   

    N_energy=photon_energy/jdos_spacing
    attenuation_layer=1.0_dp

        do atom=1,max_atoms
            attenuation_layer(N_energy,atom)=exp(-(absorp_photo(N_energy,atom)*light_path(atom))*1E-10)
            absorption_layer(N_energy,atom)=absorp_photo(N_energy,atom)*&
            thickness_atom(atom)*1E-10

        end do


    I_0=1.0_dp
    I_layer=1.0_dp
    I_layer(:,1)=I_0-reflect_photo(N_energy,1)
    if(max_layer.gt.1)then

            do atom=first_atom_second_l,max_atoms
                I_layer(N_energy,layer(atom))=I_layer(N_energy,layer(atom)-1)*&
                exp(-(absorp_photo(N_energy,atom)*light_path(atom)*1E-10))
                if(I_layer(N_energy,layer(atom)).lt.0.0_dp)then
                    I_layer(N_energy,layer(atom))=0.0_dp
                end if
             end do

    end if

    if(allocated(light_path)) then
        deallocate(light_path,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate light_path')
    end if
    if(allocated(attenuation_layer)) then
        deallocate(attenuation_layer,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate attenuation_layer')
    end if


 end subroutine calc_absorp_layer


  !***************************************************************
  subroutine calc_electron_esc
    !***************************************************************
    ! This subroutine calculates the electron escape depth

    use od_constants, only : dp,deg_to_rad
    use od_electronic, only : nbands, nspins, band_energy, efermi
    use od_cell, only : nkpoints, num_kpoints_on_node, num_atoms, &
    atoms_pos_cart_photo
    use od_io, only : io_error
    use od_comms, only : my_node_id
    use od_parameters, only : imfp_const

    integer :: atom,N,N_spin,n_eigen,ierr
    real(kind=dp) :: x,g2,efermi_scaled,a !(Evacuum-Ef)

    allocate(new_atoms_coordinates(3,max_atoms))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of new_atoms_coordinates failed')
   
    !Redefine new z coordinates where the first layer is at z=0        
    new_atoms_coordinates=atoms_pos_cart_photo
    do atom=1,max_atoms
        new_atoms_coordinates(3,atom_order(atom))=atoms_pos_cart_photo(3,atom_order(atom))-&
        (atoms_pos_cart_photo(3,atom_order(1)))
    end do

    allocate(electron_esc(nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    electron_esc=0.0_dp

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
             do atom=1,max_atoms
                if(cos(theta_arpes(n_eigen,N,N_spin)*deg_to_rad).gt.0.0_dp)then
                  electron_esc(n_eigen,N,N_spin,atom)=&
                  exp((new_atoms_coordinates(3,atom_order(atom))/&
                  cos(theta_arpes(n_eigen,N,N_spin)*deg_to_rad))/imfp_const)
                end if
             end do
          end do
       end do
    end do
  end subroutine calc_electron_esc

    !***************************************************************
  subroutine bulk_emission
    !***************************************************************
    ! This subroutine calculates the electron escape depth

    use od_constants, only : dp,deg_to_rad
    use od_electronic, only : nbands, nspins, band_energy, efermi
    use od_cell, only : nkpoints, num_kpoints_on_node, num_atoms, &
    atoms_pos_cart_photo,num_atoms,num_species
    use od_comms, only : my_node_id
    use od_parameters, only : imfp_const,photon_energy,jdos_spacing,bulk_length
    use od_jdos_utils, only : jdos_nbins, E
    use od_io, only : stdout, io_error

    real(kind=dp), dimension(:,:,:,:), allocatable :: bulk_esc_tmp
    real(kind=dp), dimension(:), allocatable :: bulk_light_tmp
     real(kind=dp), dimension(:,:,:,:), allocatable :: bulk_prob_tmp
    real(kind=dp) :: bulk_thickness
    integer :: N,N_spin,n_eigen,i,N_energy,num_layers
    integer :: atom,ierr
  
    num_layers = int((imfp_const*bulk_length)/thickness_atom(max_atoms))

    allocate(bulk_esc_tmp(nbands,num_kpoints_on_node(my_node_id),nspins,num_layers))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    bulk_esc_tmp=0.0_dp
    allocate(bulk_light_tmp(num_layers))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    bulk_light_tmp=0.0_dp
    allocate(bulk_prob_tmp(nbands,num_kpoints_on_node(my_node_id),nspins,num_layers))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    bulk_prob_tmp=0.0_dp
    allocate(bulk_prob(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_electron_esc - allocation of electron_esc failed')
    bulk_prob=0.0_dp

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
              if(cos(theta_arpes(n_eigen,N,N_spin)*deg_to_rad).gt.0.0_dp)then
                do i=1,num_layers
                    bulk_esc_tmp(n_eigen,N,N_spin,i) = &
                    (exp(((new_atoms_coordinates(3,atom_order(max_atoms))-&
                    i*thickness_atom(max_atoms)/&
                    cos(theta_arpes(n_eigen,N,N_spin)*deg_to_rad)))/imfp_const))
                end do
              end if
          end do
        end do
    end do

    bulk_light_tmp(1) = I_layer(N_energy,layer(max_atoms))*&
                exp(-(absorp_photo(N_energy,max_atoms)*thickness_atom(max_atoms)*1E-10))
    do i=2,num_layers
        bulk_light_tmp(i) = bulk_light_tmp(i-1)*& 
                exp(-(absorp_photo(N_energy,max_atoms)*i*thickness_atom(max_atoms)*1E-10))
    end do
    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
              do i=1,num_layers
                  bulk_prob_tmp(n_eigen,N,N_spin,i) = &
                  bulk_esc_tmp(n_eigen,N,N_spin,i)*bulk_light_tmp(i)
              end do
                  bulk_prob(n_eigen,N,N_spin) = &
                  sum(bulk_prob_tmp(n_eigen,N,N_spin,1:num_layers))
          end do
      end do
    end do


    if(allocated(bulk_esc_tmp)) then
        deallocate(bulk_esc_tmp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate')
    end if
    if(allocated(bulk_light_tmp)) then
        deallocate(bulk_light_tmp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate')
    end if
    if(allocated(bulk_prob_tmp)) then
        deallocate(bulk_prob_tmp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate')
    end if

    if(allocated(new_atoms_coordinates)) then
        deallocate(new_atoms_coordinates,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate new_atoms_coordinates')
    end if

   if(allocated(thickness_atom)) then
        deallocate(thickness_atom,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate thickness_atom')
    end if

  end subroutine bulk_emission

  !===============================================================================
  subroutine calc_angle
  !===============================================================================
  ! This subroutine calculates the photoemission angles theta and phi
  ! Theta: angle between the photoemitted electron and the perpendicular
  !        of the surface
  ! Phi: angle between the x and y components parallel to the surface
  ! Victor Chang, 7th February 2020
  !===============================================================================
  use od_cell, only : num_kpoints_on_node, recip_lattice, &
      cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight,num_atoms
  use od_electronic, only : nbands, nspins, band_energy, efermi,&
      electrons_per_state,band_gradient,elec_read_band_gradient,num_electrons,&
      elec_read_band_curvature, band_curvature
  use od_comms, only : my_node_id
  use od_parameters, only : work_function,photon_energy,temp, &
      elec_field,surface_area, jdos_spacing, scissor_op, &
      fixed_smearing, e_units, finite_bin_correction,adaptive_smearing, &
      hybrid_linear_grad_tol,hybrid_linear,exclude_bands, num_exclude_bands,momentum 
  use od_dos_utils, only : doslin, doslin_sub_cell_corners
  use od_algorithms, only : gaussian  
  use od_io, only : stdout, io_error, seedname, io_file_unit, stdout
  use od_jdos_utils, only : jdos_utils_calculate
  use od_jdos_utils, only : jdos_nbins, E
  use od_constants, only : hbar,ev_to_j,j_to_ev,e_mass,rad_to_deg

    integer :: N,N_spin,n_eigen,n_eigen2,atom,ierr,i,j,Gx,Gy
    integer :: angle, N_energy,transitions_den,transitions_num

    real(kind=dp),allocatable, dimension(:,:,:):: E_x
    real(kind=dp),allocatable, dimension(:,:,:):: E_y

    N_energy=photon_energy/jdos_spacing

    allocate(E_x(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    E_x = 0.0_dp

    allocate(E_y(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    E_y = 0.0_dp

    allocate(E_transverse(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    E_transverse = 0.0_dp

    allocate(E_kinetic(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    E_kinetic = 0.0_dp

    allocate(theta_arpes(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    theta_arpes = 0.0_dp

    allocate(phi_arpes(nbands,num_kpoints_on_node(my_node_id),nspins))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    phi_arpes = 0.0_dp
    
    if (index(momentum,'kp')>0) then
        call elec_read_band_gradient
        call elec_read_band_curvature
    end if
    if (index(momentum,'operator')>0) then
        call elec_read_band_gradient
    end if

    call cell_calc_kpoint_r_cart

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
             if (index(momentum,'kp')>0) then
                E_x(n_eigen,N,N_spin)=abs&
                    (0.5_dp*(1/(band_curvature(n_eigen,1,1,N,N_spin)*ev_to_j*1E-20/(hbar**2)))*&
                    (band_gradient(n_eigen,1,N,N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
                E_y(n_eigen,N,N_spin)=abs&
                    (0.5_dp*(1/(band_curvature(n_eigen,2,2,N,N_spin)*ev_to_j*1E-20/(hbar**2)))*&
                    (band_gradient(n_eigen,2,N,N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
             end if
             if (index(momentum,'crystal')>0) then
                E_x(n_eigen,N,N_spin)=(((hbar**2)/(2*e_mass))*((kpoint_r_cart(1,N)*1E+10)**2))*j_to_ev
                E_y(n_eigen,N,N_spin)=(((hbar**2)/(2*e_mass))*((kpoint_r_cart(2,N)*1E+10)**2))*j_to_ev
             end if
             if (index(momentum,'operator')>0) then
                E_x(n_eigen,N,N_spin)=abs&
                    (0.5_dp*e_mass*&
                    (band_gradient(n_eigen,1,N,N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
                E_y(n_eigen,N,N_spin)=abs&
                    (0.5_dp*e_mass*&
                    (band_gradient(n_eigen,2,N,N_spin)*(ev_to_j*1E-10/hbar))**2)*j_to_ev
             end if
             E_transverse(n_eigen,N,N_spin) = E_x(n_eigen,N,N_spin)+E_y(n_eigen,N,N_spin)
           end do
       end do
    end do

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
             phi_arpes(n_eigen,N,N_spin)=&
             atan(E_x(n_eigen,N,N_spin)/E_y(n_eigen,N,N_spin))*rad_to_deg
          end do
       end do
    end do

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
             E_kinetic(n_eigen,N,N_spin) = &
             (band_energy(n_eigen,N_spin,N)+photon_energy-evacuum_eff)

     !Calculat angle
             theta_arpes(n_eigen,N,N_spin) =&
             (acos((E_kinetic(n_eigen,N,N_spin)-E_transverse(n_eigen,N,N_spin))/&
             E_kinetic(n_eigen,N,N_spin)))*rad_to_deg
          end do
       end do
    end do


    if(allocated(E_kinetic)) then
        deallocate(E_kinetic,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

    if(allocated(kpoint_r_cart)) then
        deallocate(kpoint_r_cart,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

  end subroutine calc_angle

  !===============================================================================
  subroutine calc_three_step_model    
  !===============================================================================
  ! This subroutine calculates the QE using the thre step model.
  ! Victor Chang, 7th February 2020
  !===============================================================================

  use od_cell, only : num_kpoints_on_node, recip_lattice, &
      cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight,num_atoms
  use od_electronic, only : nbands, nspins, band_energy, efermi,&
      electrons_per_state,band_gradient,elec_read_band_gradient,num_electrons,&
      elec_read_band_curvature, band_curvature
  use od_comms, only : my_node_id
  use od_parameters, only : work_function,photon_energy, &
      elec_field,surface_area, jdos_spacing, scissor_op, &
      temp, e_units, finite_bin_correction,adaptive_smearing, &
      hybrid_linear_grad_tol,hybrid_linear,exclude_bands, num_exclude_bands,&
      write_photo_matrix
  use od_dos_utils, only : doslin, doslin_sub_cell_corners
  use od_algorithms, only : gaussian  
  use od_io, only : stdout, io_error, seedname, io_file_unit, stdout
  use od_jdos_utils, only : jdos_utils_calculate
  use od_jdos_utils, only : jdos_nbins, E
  use od_constants, only : pi, kB


    integer :: N,N_spin,n_eigen,n_eigen2,atom,ierr,i,j,Gx,Gy
    integer :: angle, N_energy
    real(kind=dp),allocatable, dimension(:,:,:,:) :: delta_temp
    real(kind=dp) :: width, norm_gaus,norm_vac,vac_g,transverse_g
    real(kind=dp) :: kbT,fermi_dirac,t_den, qe_factor,band_eff
    logical :: fixed
    integer :: matrix_unit

      width = (1.0_dp/11604.45_dp)*temp
      qe_factor=1.0_dp/(2*pi*surface_area)
      norm_vac = gaussian(0.0_dp,width,0.0_dp)

    N_energy=photon_energy/jdos_spacing

    if(allocated(epsilon)) then
       deallocate(epsilon,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon')
    end if

    if(allocated(epsilon_sum)) then
       deallocate(epsilon_sum,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon_sum')
    end if

    if(allocated(refract)) then
       deallocate(refract,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate refract')
    end if

    allocate(field_emission(nbands,nspins,num_kpoints_on_node(my_node_id)))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of enery failed')
     field_emission=0.0_dp

    if(elec_field.gt.0.0_dp)then
        call calc_field_emission
    end if

    allocate(qe_tsm(nbands,nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms+1))
    if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
    qe_tsm = 0.0_dp

    call jdos_utils_calculate_delta(delta_temp)

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
              band_eff=(band_energy(n_eigen,N_spin,N)-efermi)
              fermi_dirac=1.0_dp/(exp((band_eff/(kB*temp)))+1.0_dp)
              if((photon_energy-E_transverse(n_eigen,N,N_spin)).le.(evacuum_eff-efermi))then
                  transverse_g=gaussian((photon_energy-E_transverse(n_eigen,N,N_spin)),&
                  width,(evacuum_eff-efermi))/norm_vac
              else
                  transverse_g = 1.0_dp
              end if
              if((band_energy(n_eigen,N_spin,N)+photon_energy).lt.evacuum_eff)then
                  vac_g=gaussian((band_energy(n_eigen,N_spin,N)+photon_energy)+&
                  scissor_op,width,evacuum_eff)/norm_vac
              else
                  vac_g= 1.0_dp
              end if
                  do n_eigen2=1,nbands
                    if(band_energy(n_eigen2,N_spin,N).lt.efermi) cycle
                      do atom=1,max_atoms
                         qe_tsm(n_eigen,n_eigen2,N,N_spin,atom)=&
                         (matrix_weights(n_eigen,n_eigen2,N,N_spin,1)*&
                         delta_temp(n_eigen,n_eigen2,N,N_spin)*&
                         electron_esc(n_eigen,N,N_spin,atom)*&
                         electrons_per_state*kpoint_weight(N)*&
                         (I_layer(N_energy,layer(atom)))*&
                         qe_factor*transverse_g*vac_g*fermi_dirac*&
                         (pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin)/&
                         pdos_weights_k_band(n_eigen,N,N_spin)))*&
                         (1+field_emission(n_eigen,N_spin,N))
                      end do
                         qe_tsm(n_eigen,n_eigen2,N,N_spin,max_atoms+1)=&
                         (matrix_weights(n_eigen,n_eigen2,N,N_spin,1)*&
                         delta_temp(n_eigen,n_eigen2,N,N_spin)*&
                         bulk_prob(n_eigen,N,N_spin)*&
                         electrons_per_state*kpoint_weight(N)*&
                         qe_factor*transverse_g*vac_g*fermi_dirac*&
                         (pdos_weights_atoms(atom_order(max_atoms),n_eigen,N,N_spin)/&
                         pdos_weights_k_band(n_eigen,N,N_spin)))*&
                         (1+field_emission(n_eigen,N_spin,N))
                   end do
          end do
       end do
    end do
    if(allocated(delta_temp)) then
        deallocate(delta_temp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

    if(allocated(optical_matrix_weights)) then
        deallocate(optical_matrix_weights,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if



    if (index(write_photo_matrix,'slab')>0) then 
  call cell_calc_kpoint_r_cart


    open(unit=matrix_unit,action='write',file=trim(seedname)//'_matrix.dat')
    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
                         write(matrix_unit,*) sum(qe_tsm(n_eigen,1:nbands,N,N_spin,1:max_atoms+1)),& 
                         (kpoint_r_cart(1,N)),(kpoint_r_cart(2,N)),&
                         band_energy(n_eigen,N_spin,N)
          end do
       end do
    end do



    close(unit=matrix_unit)
    end if

  end subroutine calc_three_step_model     


  !===============================================================================
  subroutine calc_one_step_model    
  !===============================================================================
  ! This subroutine calculates the QE using a one step model.
  ! Victor Chang, 7th February 2020
  !===============================================================================

  use od_cell, only : num_kpoints_on_node, recip_lattice,&
      cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight,num_atoms
  use od_electronic, only : nbands, nspins, band_energy, efermi,&
      electrons_per_state,band_gradient,elec_read_band_gradient,num_electrons,&
      elec_read_band_curvature, band_curvature
  use od_comms, only : my_node_id
  use od_parameters, only : work_function,photon_energy, &
      elec_field,surface_area, jdos_spacing, scissor_op, &
      e_units, finite_bin_correction,adaptive_smearing, &
      hybrid_linear_grad_tol,hybrid_linear,exclude_bands, num_exclude_bands,&
      temp,write_photo_matrix
  use od_dos_utils, only : doslin, doslin_sub_cell_corners
  use od_algorithms, only : gaussian  
  use od_io, only : stdout, io_error, seedname, io_file_unit, stdout
  use od_jdos_utils, only : jdos_utils_calculate
  use od_jdos_utils, only : jdos_nbins, E
  use od_constants, only : pi, kB


    integer :: N,N_spin,n_eigen,n_eigen2,atom,ierr,i,j,Gx,Gy
    integer :: angle, N_energy,transitions_den,transitions_num
    real(kind=dp),allocatable, dimension(:,:,:,:) :: delta_temp
    integer :: matrix_unit
    real(kind=dp) :: width, norm_gaus,norm_vac,vac_g,transverse_g
    real(kind=dp) :: fermi_dirac,t_den, qe_factor,band_eff
    logical :: fixed

     width = (1.0_dp/11604.45_dp)*temp
     qe_factor=1.0_dp/(2*pi*surface_area)
     norm_vac = gaussian(0.0_dp,width,0.0_dp)

     N_energy=photon_energy/jdos_spacing

     if(allocated(epsilon)) then
        deallocate(epsilon,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon')
     end if

     if(allocated(epsilon_sum)) then
        deallocate(epsilon_sum,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon_sum')
     end if

     if(allocated(refract)) then
        deallocate(refract,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate refract')
     end if

     if(allocated(matrix_weights)) then
        deallocate(matrix_weights,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate matrix_weights')
     end if

    allocate(field_emission(nbands,nspins,num_kpoints_on_node(my_node_id)))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of enery failed')
     field_emission=0.0_dp

      if(elec_field.gt.0.0_dp)then
          call calc_field_emission
      end if

     allocate(qe_osm(nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms+1))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
     qe_osm = 0.0_dp

     do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
              band_eff=(band_energy(n_eigen,N_spin,N)-efermi)
              fermi_dirac=1.0_dp/(exp((band_eff/(kB*temp)))+1.0_dp)
            if((photon_energy-E_transverse(n_eigen,N,N_spin)).le.(evacuum_eff-efermi))then
              transverse_g=gaussian((photon_energy-E_transverse(n_eigen,N,N_spin)),&
              width,(evacuum_eff-efermi))/norm_vac
            else
              transverse_g = 1.0_dp
            end if
            if((band_energy(n_eigen,N_spin,N)+photon_energy).lt.evacuum_eff)then
              vac_g=gaussian((band_energy(n_eigen,N_spin,N)+photon_energy)+&
              scissor_op,width,evacuum_eff)/norm_vac
            else
              vac_g= 1.0_dp
            end if
                    n_eigen2=nbands+1
                    do atom=1,max_atoms
                      qe_osm(n_eigen,N,N_spin,atom)=&
                      (foptical_matrix_weights(n_eigen,n_eigen2,N,N_spin,1)*&
                      (electron_esc(n_eigen,N,N_spin,atom))*&
                      electrons_per_state*kpoint_weight(N)*&
                      (I_layer(N_energy,layer(atom)))*&
                      qe_factor*&
                      transverse_g*vac_g*fermi_dirac*&
                      (pdos_weights_atoms(atom_order(atom),n_eigen,N,N_spin)/&
                      pdos_weights_k_band(n_eigen,N,N_spin)))*&
                      (1+field_emission(n_eigen,N_spin,N))
                   end do
                   qe_osm(n_eigen,N,N_spin,max_atoms+1)=&
                      (foptical_matrix_weights(n_eigen,n_eigen2,N,N_spin,1)*&
                      bulk_prob(n_eigen,N,N_spin)*&
                      electrons_per_state*kpoint_weight(N)*&
                      qe_factor*transverse_g*vac_g*fermi_dirac*&
                      (pdos_weights_atoms(atom_order(max_atoms),n_eigen,N,N_spin)/&
                      pdos_weights_k_band(n_eigen,N,N_spin)))*&!+&
                      (1+field_emission(n_eigen,N_spin,N))
          end do
       end do
    end do

   if (index(write_photo_matrix,'slab')>0) then
   call cell_calc_kpoint_r_cart

    open(unit=matrix_unit,action='write',file=trim(seedname)//'_matrix.dat')
    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
                         write(matrix_unit,*)sum(qe_osm(n_eigen,N,N_spin,1:max_atoms)),&
                         (kpoint_r_cart(1,N)),(kpoint_r_cart(2,N)),&
                         band_energy(n_eigen,N_spin,N)
          end do
       end do
    end do

    close(unit=matrix_unit)
    end if

    if(allocated(foptical_matrix_weights)) then
        deallocate(foptical_matrix_weights,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

end subroutine calc_one_step_model     

  !===============================================================================
  subroutine jdos_utils_calculate_delta(delta_temp)
  !===============================================================================
  ! It is required to evaluate the delta funcion.
  ! Victor Chang, 7th February 2020
  !===============================================================================
    use od_parameters, only : linear, fixed, adaptive, quad, iprint, dos_per_volume
    use od_electronic, only : elec_read_band_gradient,band_gradient,nspins,electrons_per_state,&
         num_electrons,efermi_set
    use od_comms, only      : on_root
    use od_io, only         : stdout,io_error,io_time
    use od_cell, only       : cell_volume
    use od_dos_utils, only     : dos_utils_set_efermi
    use od_jdos_utils, only : setup_energy_scale

    implicit none
    integer :: ierr
    real(kind=dp) :: time0, time1

    real(kind=dp),intent(out), allocatable, optional    :: delta_temp(:,:,:,:)  !I've added this
    real(kind=dp), allocatable :: jdos_adaptive(:,:)
    real(kind=dp), allocatable :: jdos_fixed(:,:)
    real(kind=dp), allocatable :: jdos_linear(:,:)

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S
    ! If we're using one of the more accurate roadening schemes we also need to read in the
    ! band gradients too
    if(quad.or.linear.or.adaptive) then
       if(.not.allocated(band_gradient)) call elec_read_band_gradient
    endif
    !-------------------------------------------------------------------------------

    if(.not.efermi_set) call dos_utils_set_efermi

    time0=io_time()

    call setup_energy_scale(E)

    if(fixed)then
              call calculate_delta('f',delta_temp)
    endif
    if(adaptive)then
              call calculate_delta('a',delta_temp)

    endif
    if(linear)then
          call calculate_delta('l',delta_temp)
    endif

    if(quad) then
       call io_error("quadratic broadening not implemented")
    endif

    time1=io_time()
    if(on_root.and.iprint>1) then
       write(stdout,'(1x,a59,f11.3,a8)') &
            '+ Time to calculate Joint Density of States              &
            &      ',time1-time0,' (sec) +'
    endif
    !-------------------------------------------------------------------------------

       if(dos_per_volume) then
       if(fixed) then
          jdos_fixed=jdos_fixed/cell_volume
       endif
       if(adaptive) then
          jdos_adaptive=jdos_adaptive/cell_volume
       endif
       if(linear) then
          jdos_linear=jdos_linear/cell_volume
       endif

       ! if(quad) then
       !    dos_quad=dos_quad/cell_volume
       !    intdos_quad=intdos_quad/cell_volume
       ! endif
    endif

  end subroutine jdos_utils_calculate_delta

    !===============================================================================
  subroutine calculate_delta(delta_type, delta_temp)
    !===============================================================================
    ! This subroutine evaluates the delta function between the valence band
    ! and the conduction band using the method specified in the input.
    ! Victor Chang, 7 February 2020
    !===============================================================================
    use od_comms, only : my_node_id, on_root
    use od_cell, only : num_kpoints_on_node, kpoint_grid_dim,kpoint_weight,&
         &recip_lattice, num_atoms
    use od_parameters, only : adaptive_smearing, fixed_smearing, iprint, &
         &finite_bin_correction, scissor_op,hybrid_linear_grad_tol,hybrid_linear, exclude_bands, num_exclude_bands,&
         photo,photon_energy,jdos_spacing
    use od_io, only : io_error,stdout
    use od_electronic, only : band_gradient,nbands,band_energy,nspins,electrons_per_state, &
         & efermi
    use od_dos_utils, only : doslin, doslin_sub_cell_corners
    use od_algorithms, only : gaussian
    implicit none

    integer :: ik,is,ib,idos,jb,i
    integer :: N,N_spin,n_eigen,n_eigen2,atom
    integer :: N2,N_geom, ierr
    real(kind=dp) :: cuml, width, adaptive_smearing_temp, dos_test
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4), sub_cell_length(1:3)

    character(len=1), intent(in)                      :: delta_type
    real(kind=dp),intent(inout),allocatable, optional :: delta_temp(:,:,:,:)
   logical :: linear,fixed,adaptive,force_adaptive

    linear=.false.
    fixed=.false.
    adaptive=.false.

    select case (delta_type)
    case ("l")
       linear=.true.
    case("a")
       adaptive=.true.
    case("f")
       fixed=.true.
    case default
       call io_error (" ERROR : unknown jdos_type in jcalculate_dos ")
    end select

    width=0.0_dp


    if(linear.or.adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:),dp)/2.0_dp
    if(adaptive.or.hybrid_linear) then
       do i= 1,3
          sub_cell_length(i)=sqrt(recip_lattice(i,1)**2+recip_lattice(i,2)**2+recip_lattice(i,3)**2)*step(i)
       enddo
       adaptive_smearing_temp=adaptive_smearing*sum(sub_cell_length)/3.0_dp
    endif

    if(fixed) width=fixed_smearing




       allocate(delta_temp(nbands,nbands,num_kpoints_on_node(my_node_id), nspins),stat=ierr)
       if(ierr/=0) call io_error('Error: calculate_jdos - failed to allocate weighted_jdos')
       delta_temp=0.0_dp
    if(iprint>1 .and. on_root) then
       write(stdout,'(1x,a78)') '+------------------------------ Calculate JDOS ------------------------------+'
    endif

    do ik=1,num_kpoints_on_node(my_node_id)
       if(iprint>1 .and. on_root) then
          if (mod(real(ik,dp),10.0_dp) == 0.0_dp) write(stdout,'(1x,a1,a38,i4,a3,i4,1x,a14,3x,a10)') ',', &
               &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id),'on this node.',"<-- JDOS |"
       endif
       do is=1,nspins
          occ_states: do ib=1,nbands
             if(num_exclude_bands>0) then
                if(any(exclude_bands==ib)) cycle
             endif
             if(band_energy(ib,is,ik).ge.efermi) cycle occ_states
             unocc_states : do jb=1,nbands
                if(band_energy(jb,is,ik).lt.efermi) cycle unocc_states
                if(linear.or.adaptive) grad(:) = band_gradient(jb,:,ik,is)-band_gradient(ib,:,ik,is)

                ! If the band is very flat linear broadening can have problems describing it. In this case, fall back to
                ! adaptive smearing (and take advantage of FBCS if required).
                force_adaptive=.false.
                if(hybrid_linear.and.(hybrid_linear_grad_tol>sqrt(dot_product(grad,grad)))) force_adaptive=.true.
                if(linear.and..not.force_adaptive) call doslin_sub_cell_corners(grad,step,band_energy(jb,is,ik)-&
                                                        &band_energy(ib,is,ik)+scissor_op,EV)
                if(adaptive.or.force_adaptive) width = sqrt(dot_product(grad,grad))*adaptive_smearing_temp

                ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
                ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
                if(finite_bin_correction.and.(width<delta_bins)) width = delta_bins

                   idos = photon_energy/jdos_spacing
                   ! The linear method has a special way to calculate the integrated dos
                   ! we have to take account for this here.
                    if(linear.and..not.force_adaptive)then
                      delta_temp(ib,jb,ik,is)=doslin(EV(0),EV(1),EV(2),EV(3),EV(4),E(idos),cuml)

                    else
                       delta_temp(ib,jb,ik,is)=gaussian((band_energy(jb,is,ik)-band_energy(ib,is,ik))+scissor_op,width,E(idos))!&
                    endif

             end do unocc_states
          end do occ_states
       end do
    end do


    if(iprint>1 .and. on_root) then
       write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    endif

  end subroutine calculate_delta

    !===============================================================================
    subroutine weighted_mean_te
    !===============================================================================    
    ! This subroutine calculates the weighted arithmetic mean transverse energy
    ! sum(QE*mte)/(total QE)
    ! Victor Chang, 7 February 2020
    !===============================================================================  
  use od_cell, only : num_kpoints_on_node, recip_lattice, &
      cell_calc_kpoint_r_cart, kpoint_r_cart, kpoint_weight,num_atoms,&
      num_atoms, atoms_pos_cart_photo, atoms_label_tmp
  use od_electronic, only : nbands, nspins, band_energy, efermi,&
      electrons_per_state,band_gradient,elec_read_band_gradient,num_electrons,&
      elec_read_band_curvature, band_curvature
  use od_comms, only : my_node_id
  use od_parameters, only : work_function,photon_energy, &
      elec_field, photo_model,exclude_bands, num_exclude_bands
  use od_dos_utils, only : doslin, doslin_sub_cell_corners
  use od_algorithms, only : gaussian  
  use od_io, only : stdout, io_error, seedname, io_file_unit, stdout
  use od_jdos_utils, only : jdos_utils_calculate
  use od_jdos_utils, only : jdos_nbins, E

    integer :: N,N_spin,n_eigen,n_eigen2,atom,ierr,i,j,Gx,Gy
    integer :: angle, N_energy,transitions_den,transitions_num
    real(kind=dp) :: mean_te
    real(kind=dp),allocatable, dimension(:,:,:,:,:) :: te_tsm_temp
    real(kind=dp),allocatable, dimension(:,:,:,:) :: te_osm_temp
    real(kind=dp) :: width, norm_gaus,norm_vac,vac_g,transverse_g
    real(kind=dp) :: kbT,fermi_dirac,t_den, qe_factor
    logical :: fixed

     if (index(photo_model,'3step')>0) then
       allocate(te_tsm_temp(nbands,nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms+1))
       if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
       te_tsm_temp = 0.0_dp

       do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin=1,nspins                    ! Loop over spins
             do n_eigen=1,nbands
                do n_eigen2=1,nbands
                   if(band_energy(n_eigen2,N_spin,N).lt.efermi) cycle
                   do atom=1,max_atoms
                      te_tsm_temp(n_eigen,n_eigen2,N,N_spin,atom)=&
                      E_transverse(n_eigen,N,N_spin)*qe_tsm(n_eigen,n_eigen2,N,N_spin,atom)
                   end do
                     te_tsm_temp(n_eigen,n_eigen2,N,N_spin,max_atoms+1)=&
                     E_transverse(n_eigen,N,N_spin)*qe_tsm(n_eigen,n_eigen2,N,N_spin,max_atoms+1)
                end do
             end do
          end do
       end do

       mean_te=sum(te_tsm_temp(1:nbands,1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))/&
       sum(qe_tsm(1:nbands,1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))

       if(allocated(te_tsm_temp)) then
          deallocate(te_tsm_temp,stat=ierr)
          if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
       end if
     
    end if

    if (index(photo_model,'1step')>0) then

      allocate(te_osm_temp(nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms+1))
      if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_numerator failed')
      te_osm_temp = 0.0_dp
 
       do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin=1,nspins                    ! Loop over spins
             do n_eigen=1,nbands
!                if(band_energy(n_eigen,N_spin,N).ge.efermi) cycle
                do atom=1,max_atoms
                   te_osm_temp(n_eigen,N,N_spin,atom)=&
                   E_transverse(n_eigen,N,N_spin)*qe_osm(n_eigen,N,N_spin,atom)
                end do
                   te_osm_temp(n_eigen,N,N_spin,max_atoms+1)=&
                   E_transverse(n_eigen,N,N_spin)*qe_osm(n_eigen,N,N_spin,max_atoms+1)
             end do
          end do
       end do

       mean_te=sum(te_osm_temp(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))/&
       sum(qe_osm(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))

       if(allocated(te_osm_temp)) then
          deallocate(te_osm_temp,stat=ierr)
          if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
       end if

    end if
    if (index(photo_model,'3step')>0) then
    write(stdout,'(1x,a78)') '+---------------------------- Photoemission ---------------------------------+'
    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    write(stdout,'(1x,a15,f15.4,1x,a25,f15.4,1x,a6)')'| Work Function',work_function,&
                                                     'eV          Photon Energy',photon_energy,'eV   |'
    write(stdout,'(a26,f10.4,1x,a26,f10.4,a6)')'| Effective Work Function',work_function_eff,&
                                                     ' eV        Electric field',elec_field,'V/A  |'
    write(stdout,'(1x,a78)') '| Final state : Bloch state                                                  |'
    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    write(stdout,'(1x,a78)') '| Atom |  Atom Order  |   Layer   |             Quantum Efficiency           |'

    do atom=1,max_atoms
        write(stdout,*)"|  ",trim(atoms_label_tmp(atom_order(atom))),atom_order(atom),&
                   layer(atom),'              ',&
                   sum(qe_tsm(1:nbands,1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,atom)),"      |"
    end do
       write(stdout,*)"|  Bulk                                     ",&
             sum(qe_tsm(1:nbands,1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,max_atoms+1)),"      |"

       write(stdout,*) '|  Total quantum efficiency (electrons/photon):',&
       sum(qe_tsm(1:nbands,1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1)),'   |'
    end if

    if (index(photo_model,'1step')>0) then
    write(stdout,'(1x,a78)') '+---------------------------- Photoemission ---------------------------------+'
    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    write(stdout,'(1x,a15,f15.4,1x,a25,f15.4,1x,a5)')'| Work Function',work_function,&
                                                     'eV          PhotonEnergy',photon_energy,'eV  |'
    write(stdout,'(a26,f10.4,1x,a26,f10.4,a6)')'| Effective Work Function',work_function_eff,&
                                                     ' eV        Electric field',elec_field,'V/A  |'
    write(stdout,'(1x,a78)') '| Final state : Free electron state                                          |'
    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    write(stdout,'(1x,a78)') '| Atom |  Atom Order  |   Layer   |             Quantum Efficiency           |'

    do atom=1,max_atoms
      write(stdout,*)"|  ",trim(atoms_label_tmp(atom_order(atom))),atom_order(atom),&
                   layer(atom),'              ',&
                   sum(qe_osm(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,atom)),"      |"
    end do
      write(stdout,*)"|  Bulk                                     ",&
                   sum(qe_osm(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,max_atoms+1)),"      |"

      write(stdout,*)'|  Total quantum efficiency (electrons/photon):',&
                      sum(qe_osm(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1)),'   |'

    end if
    write(stdout,*) '|  Weighted mean transverse energy (eV):',mean_te,'          |'
    if(elec_field.gt.0.0_dp)then
    write(stdout,*) '|  Total field emission (electrons/A^2):',total_field_emission, '          |'
    end if
    write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'



    end subroutine weighted_mean_te

  subroutine binding_energy_spread
  ! This subroutine applies a Gaussian broadenning to the binding energy
  ! Additionally, it takes the photoemission angles theta and phi as inputs
  ! Victor Chang, 7 February 2020

    use od_cell, only : num_kpoints_on_node,cell_calc_kpoint_r_cart,kpoint_r_cart
    use od_electronic, only : nbands, nspins,band_energy, efermi
    use od_parameters, only : work_function,photon_energy,fixed_smearing,&
    photo_model,theta_lower,theta_upper,phi_lower,phi_upper
    use od_algorithms, only : gaussian
    use od_comms, only : my_node_id
   use od_io, only : stdout, io_error, seedname, io_file_unit

    real(kind=dp) :: norm_trans,qe_norm
    real(kind=dp), allocatable, dimension(:,:,:,:) :: binding_temp
     real(kind=dp), allocatable, dimension(:,:,:,:) :: qe_temp

    integer :: e_scale,ierr
    integer :: N,N_spin,n_eigen,n_eigen2,atom

    integer :: test_unit

    max_energy=int((photon_energy-work_function)*1000)+100





     if(allocated(E_transverse)) then
        deallocate(E_transverse,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
     end if

    allocate(t_energy(max_energy))
    t_energy = 0.0_dp

    allocate(weighted_temp(max_energy,num_kpoints_on_node(my_node_id),nspins,nbands,max_atoms+1))
    weighted_temp = 0.0_dp

    allocate(qe_temp(nbands,num_kpoints_on_node(my_node_id),nspins,max_atoms+1))
    qe_temp = 0.0_dp

    allocate(binding_temp(max_energy,num_kpoints_on_node(my_node_id),nspins,nbands))
    binding_temp=0.0_dp

    do e_scale=1,max_energy
       t_energy(e_scale) = real(e_scale-1,dp)/1000
    end do

    do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
       do N_spin=1,nspins                    ! Loop over spins
          do n_eigen=1,nbands
             do e_scale=1,max_energy
                binding_temp(e_scale,N,N_spin,n_eigen)=&
                gaussian((efermi-band_energy(n_eigen,N_spin,N)),fixed_smearing,t_energy(e_scale)) 
             end do
          end do
       end do
    end do


  if (index(photo_model,'3step')>0) then
    do atom=1,max_atoms+1
       do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
          do N_spin=1,nspins                    ! Loop over spins
             do n_eigen=1,nbands
                qe_temp(n_eigen,N,N_spin,atom)=sum(qe_tsm(n_eigen,1:nbands,N,N_spin,atom))
             end do
          end do
       end do
    end do

    if(allocated(qe_tsm)) then
        deallocate(qe_tsm,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate qe matrix')
    end if

    do e_scale=1,max_energy
       do atom=1,max_atoms+1
          do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
             do N_spin=1,nspins                    ! Loop over spins
                do n_eigen=1,nbands
                   if(theta_arpes(n_eigen,N,N_spin).ge.theta_lower.and.&
                     theta_arpes(n_eigen,N,N_spin).le.theta_upper)then
                       if(phi_arpes(n_eigen,N,N_spin).ge.phi_lower.and.&
                         phi_arpes(n_eigen,N,N_spin).le.phi_upper)then
                           weighted_temp(e_scale,N,N_spin,n_eigen,atom)=&
                           binding_temp(e_scale,N,N_spin,n_eigen)*qe_temp(n_eigen,N,N_spin,atom)
                       end if
                    end if
                 end do
              end do
           end do
        end do
     end do

     qe_norm = sum(qe_temp(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))&
     /sum(weighted_temp(1:max_energy,1:num_kpoints_on_node(my_node_id),1:nspins,1:nbands,1:max_atoms+1))

     do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
              do e_scale=1,max_energy
                 do atom=1,max_atoms+1
                    weighted_temp(e_scale,N,N_spin,n_eigen,atom)=&
                    weighted_temp(e_scale,N,N_spin,n_eigen,atom)*qe_norm
                 end do
              end do
           end do
        end do
     end do

     if(allocated(qe_temp)) then
        deallocate(qe_temp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
     end if
   end if
   
   if (index(photo_model,'1step')>0) then
     do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
              if(theta_arpes(n_eigen,N,N_spin).ge.theta_lower.and.&
                theta_arpes(n_eigen,N,N_spin).le.theta_upper)then
                if(phi_arpes(n_eigen,N,N_spin).ge.phi_lower.and.&
                  phi_arpes(n_eigen,N,N_spin).le.phi_upper)then
!                  if(band_energy(n_eigen,N_spin,N).ge.efermi) cycle
                     do e_scale=1,max_energy
                        do atom=1,max_atoms+1
                           weighted_temp(e_scale,N,N_spin,n_eigen,atom)=&
                           binding_temp(e_scale,N,N_spin,n_eigen)*qe_osm(n_eigen,N,N_spin,atom) 
                        end do
                     end do
                  end if
               end if
            end do
         end do
      end do

      qe_norm = sum(qe_osm(1:nbands,1:num_kpoints_on_node(my_node_id),1:nspins,1:max_atoms+1))&
         /sum(weighted_temp(1:max_energy,1:num_kpoints_on_node(my_node_id),1:nspins,1:nbands,1:max_atoms+1))

      do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
         do N_spin=1,nspins                    ! Loop over spins
            do n_eigen=1,nbands
               do e_scale=1,max_energy
                  do atom=1,max_atoms+1
                     weighted_temp(e_scale,N,N_spin,n_eigen,atom)=&
                     weighted_temp(e_scale,N,N_spin,n_eigen,atom)*qe_norm
                  end do
               end do
            end do
         end do
      end do

      if(allocated(qe_osm)) then
        deallocate(qe_osm,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
      end if
    end if

     if(allocated(binding_temp)) then
        deallocate(binding_temp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
     end if



   end subroutine binding_energy_spread

    !***************************************************************
   subroutine write_qe_output_files
    !***************************************************************
    ! This subroutine write either the transverse energy or the binding energy
    ! after the Gaussian broadening has been applied.
    ! Victor Chang, 7 February 2020

   use od_cell, only : num_kpoints_on_node
   use od_electronic, only : nbands, nspins,efermi,band_energy
   use od_comms, only : my_node_id
   use od_io, only : stdout, io_error, seedname, io_file_unit
   use od_parameters, only : photo_model


    integer :: transverse_unit, binding_unit
    integer :: N,N_spin,n_eigen,n_eigen2,atom,ierr,e_scale

    real(kind=dp),allocatable, dimension(:,:) :: qe_atom


     allocate(qe_atom(max_energy,max_atoms+1))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of qe_atom failed')
     qe_atom=0.0_dp

     do e_scale=1,max_energy !loop over transverse energy
        do atom=1,max_atoms+1
           qe_atom(e_scale,atom) = &
           sum(weighted_temp(e_scale,1:num_kpoints_on_node(my_node_id),1:nspins,1:nbands,atom))
        end do
     end do

    open(unit=binding_unit,action='write',file=trim(seedname)//'_binding_energy.dat')
       do e_scale=1,max_energy
          write(binding_unit,*)t_energy(e_scale),sum(qe_atom(e_scale,1:max_atoms+1)),&
                               qe_atom(e_scale,1:max_atoms+1)
       end do
    close(unit=binding_unit)


    if(allocated(weighted_temp)) then
        deallocate(weighted_temp,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

        if(allocated(qe_atom)) then
        deallocate(qe_atom,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

        if(allocated(t_energy)) then
        deallocate(t_energy,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

   end subroutine write_qe_output_files

     !***************************************************************
   subroutine effect_wf
     !***************************************************************
     !elec_field given in eV/A

     use od_parameters, only : work_function, elec_field
     use od_electronic, only : efermi
     use od_constants, only : pi, epsilon_zero


    real(kind=dp) :: z

    z=sqrt((1/(16*pi*epsilon_zero*1E-4))/elec_field)

    work_function_eff = work_function - elec_field*z -(1/(16*pi*epsilon_zero*1E-4))/z

    evacuum_eff = work_function_eff + efermi



   end subroutine effect_wf

     !***************************************************************
   subroutine calc_field_emission
     !***************************************************************
     ! This subroutine calculates the Schottky effect
     !elec_field given in V/m

     use od_cell, only : num_kpoints_on_node
     use od_parameters, only : work_function, elec_field,temp,surface_area
     use od_electronic, only : efermi,band_energy,nbands, nspins
     use od_io, only : stdout, io_error
     use od_comms, only : my_node_id
     use od_constants, only : pi,epsilon_zero,ge,kB

    integer :: ierr
    real(kind=dp),allocatable, dimension(:,:,:) :: field_energy
    real(kind=dp),allocatable, dimension(:,:,:) :: tunnel_prob
    real(kind=dp),allocatable, dimension(:,:,:) :: G
    real(kind=dp),allocatable, dimension(:,:,:) :: temp_emission
    real(kind=dp) :: fermi_dirac,band_eff

    integer :: N,N_spin,n_eigen,z_distance,z,z_max
    real(kind=dp) :: lower_zero,upper_zero,z_range

     evacuum = efermi+work_function
     allocate(field_energy(nbands,nspins,num_kpoints_on_node(my_node_id)))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of enery failed')
     field_energy=0.0_dp
     allocate(G(nbands,nspins,num_kpoints_on_node(my_node_id)))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of enery failed')
     G=0.0_dp
     allocate(temp_emission(nbands,nspins,num_kpoints_on_node(my_node_id)))
     if(ierr/=0) call io_error('Error: calc_quantum_efficiency - allocation of enery failed')
     temp_emission=0.0_dp

     do N=1,num_kpoints_on_node(my_node_id)   ! Loop over kpoints
        do N_spin=1,nspins                    ! Loop over spins
           do n_eigen=1,nbands
              band_eff=(band_energy(n_eigen,N_spin,N)-efermi)
              fermi_dirac=1.0_dp/(exp((band_eff/(kB*temp)))+1.0_dp)
              field_energy(n_eigen,N_spin,N) = abs(evacuum - band_energy(n_eigen,N_spin,N)) 

              if((elec_field/(4*pi*epsilon_zero*1E-4)*elec_field).lt.(field_energy(n_eigen,N_spin,N)**2))then

                   !find the zeros
                   lower_zero=(field_energy(n_eigen,N_spin,N)-&
                   sqrt(field_energy(n_eigen,N_spin,N)**2-&
                   (elec_field/(4*pi*epsilon_zero*1E-4))))/&
                   (2*elec_field)

                   upper_zero=(field_energy(n_eigen,N_spin,N)+&
                   sqrt(field_energy(n_eigen,N_spin,N)**2-&
                   (elec_field/(4*pi*epsilon_zero*1E-4))))/&
                   (2*elec_field)

                   z_range = upper_zero-lower_zero
                 
                   if(z_range.eq.z_range)then
                       z_max = int(z_range)/0.5
                       do z=1,z_max-1
                           G(n_eigen,N_spin,N) = G(n_eigen,N_spin,N) + &
                               (sqrt(field_energy(n_eigen,N_spin,N)-&
                               elec_field*((z*0.5)+lower_zero)-&
                               (1/(16*pi*epsilon_zero*1E-4)/((z*0.5)+lower_zero))))
                       end do
                    end if
                    field_emission(n_eigen,N_spin,N) = 1/(1+exp(G(n_eigen,N_spin,N)))
              end if
               temp_emission(n_eigen,N_spin,N) = field_emission(n_eigen,N_spin,N)*fermi_dirac

              if(temp_emission(n_eigen,N_spin,N).gt.0.0_dp)then
          end if
        end do
      end do
    end do

    total_field_emission = sum(temp_emission(1:nbands,1:nspins,1:num_kpoints_on_node(my_node_id)))/&
                       surface_area

    if(allocated(field_energy)) then
        deallocate(field_energy,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate energy')
    end if

        if(allocated(G)) then
        deallocate(G,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate G')
    end if

        if(allocated(temp_emission)) then
        deallocate(temp_emission,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate G')
    end if

     !ADD COMMENT 
   end subroutine calc_field_emission



    !***************************************************************
  subroutine photo_deallocate        
    !***************************************************************
    ! This subroutine deallocates all the quantities which have not 
    ! been deallocated yet
    
     use od_io, only : stdout, io_error
    integer :: ierr


    if(allocated(pdos_weights_atoms)) then
       deallocate(pdos_weights_atoms,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_atoms')
    end if

    if(allocated(pdos_weights_k_band)) then
       deallocate(pdos_weights_k_band,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate pdos_weights_k_band')
    end if

    if(allocated(optical_matrix_weights)) then
       deallocate(optical_matrix_weights,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate optical_matrix_weights')
    end if

    if(allocated(matrix_weights)) then
       deallocate(matrix_weights,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate matrix_weights')
    end if

    if(allocated(epsilon)) then
       deallocate(epsilon,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon')
    end if

    if(allocated(epsilon_sum)) then
       deallocate(epsilon_sum,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate epsilon_sum')
    end if

    if(allocated(refract)) then
       deallocate(refract,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate refract')
    end if


    if(allocated(absorp)) then
       deallocate(absorp,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate absorp')
    end if

    if(allocated(electron_esc)) then
       deallocate(electron_esc,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate electron_esc')
    end if

    if(allocated(layer)) then
       deallocate(layer,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate layer')
    end if

    if(allocated(imfp_val)) then
       deallocate(imfp_val,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate imfp_val')
    end if

    if(allocated(reflect)) then
       deallocate(reflect,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate reflect')
    end if

    if(allocated(total_transmittance)) then
       deallocate(total_transmittance,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate total_transmittance')
    end if

    if(allocated(total_absorption)) then
       deallocate(total_absorption,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate total_absorption')
    end if

    if(allocated(atom_order)) then
       deallocate(atom_order,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate atom_order')
    end if

    if(allocated(atoms_per_layer)) then
       deallocate(atoms_per_layer,stat=ierr)
       if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate atoms_per_layer')
    end if

    if(allocated(E_transverse)) then
        deallocate(E_transverse,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if
    if(allocated(phi_arpes)) then
        deallocate(phi_arpes,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if

    if(allocated(theta_arpes)) then
        deallocate(theta_arpes,stat=ierr)
        if(ierr/=0) call io_error('Error: photo_deallocate - failed to deallocate kpoint_r_cart')
    end if


  end subroutine photo_deallocate        

endmodule od_photo
