!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=========================================================================!
! E L E C T R O N I C 
! Stores variables to do with the electrons and energy eigenvalues in the
! system. At this stage I don't see it contining many functions, but it
! was important for the dos module not to have too many global variables
! AJM Dec 2010
!=========================================================================!

module od_electronic
  use od_constants, only : dp

  implicit none

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S 
  real(kind=dp), allocatable, public, save  :: band_energy(:,:,:)
  complex(kind=dp), allocatable, public, save  :: band_gradient(:,:,:,:,:)  !I've changed this from real to complex
  complex(kind=dp), allocatable, public, save  :: elnes_mat(:,:,:,:,:)

  real(kind=dp), public, save :: efermi ! The fermi energy we finally decide on
  real(kind=dp), public, save :: efermi_castep ! Fermi energy as reported by CASTEP

  real(kind=dp), allocatable, public, save :: num_electrons(:) ! Holds up-spin and 
  ! down-spin

  integer, public, save       :: nbands, nspins
  real(kind=dp), public, save :: electrons_per_state ! 2 for non-spin-P
  ! 1 for spin-P
  logical, public, save       :: spin_polarised



  type, public ::  matrix_weights_array_boundaries
     integer :: norbitals
     integer :: nbands
     integer :: nkpoints
     integer :: nspins
  end type matrix_weights_array_boundaries

  type, public :: orbitals
     integer          :: ion_no          ! Unique ion number
     integer          :: species_no      ! Unique species number
     integer          :: rank_in_species ! Unique ion number within species
     integer          :: am_channel      ! The angular momentum Channel (l)
     integer          :: shell           ! Principal quantum number (n) !n.b typically only know this for core states
     logical          :: calc_pdos       ! Should pDoS be calculated for this orbital?                                         
     character(len=1) :: am_channel_name ! Name of angular momentum channel s,p,d, etc
  end type orbitals

  type(orbitals), public, allocatable, save :: pdos_orbital(:)
  real(kind=dp), public, allocatable, save  :: pdos_weights(:,:,:,:)
  type(matrix_weights_array_boundaries), public, save :: pdos_mwab

  type(orbitals), public, allocatable, save :: elnes_orbital(:)
  type(matrix_weights_array_boundaries), public, save :: elnes_mwab


  !-------------------------------------------------------------------------!

  private

  !-------------------------------------------------------------------------!
  ! G L O B A L   F U N C T I O N S
  public :: elec_report_parameters
  public :: elec_read_band_energy
  public :: elec_read_band_gradient
  public :: elec_read_elnes_mat
  public :: elec_pdos_read
  public :: elec_pdos_read_orbitals
  public :: elec_dealloc_elnes
  public :: elec_dealloc_pdos
  public :: elec_dealloc_band_gradient

  !-------------------------------------------------------------------------! 

contains

  !=========================================================================
  subroutine elec_report_parameters
    !=========================================================================
    ! Report the electronic properties in the calculation
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables nbands,num_electrons,nkpoints,kpoint_grid_dim                                      
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_io,   only : stdout 
    use od_cell, only : kpoint_grid_dim, nkpoints
    implicit none
    integer :: i

    write (stdout,*)

    write(stdout,'(1x,a78)')    '+----------------------- Electronic Data ------------------------------------+'

    write(stdout,'(1x,a46,i14,a18)')   '|  Number of Bands                           :',nbands,"|"
    write(stdout,'(1x,a46,6x,i3,1x,a1,i3,1x,a1,i3,12x,a1)') '|  Grid size                                 :'&
         ,kpoint_grid_dim(1),'x',kpoint_grid_dim(2),'x',kpoint_grid_dim(3),'|'
    write(stdout,'(1x,a46,i14,a18)')   '|  Number of K-points                        :',nkpoints,"|"

    if (nspins > 1) then
       write(stdout,'(1x,a78)')   '|  Spin-Polarised Calculation                :           True                |'
       write(stdout,'(1x,a46,f17.2,a15)')"|  Number of up-spin electrons               :", num_electrons(1),"|"
       write(stdout,'(1x,a46,f17.2,a15)')"|  Number of down-spin electrons             :", num_electrons(2),"|"
    else
       write(stdout,'(1x,a78)')          '|  Spin-Polarised Calculation                :           False               |'
       write(stdout,'(1x,a46,f17.2,a15)')"|  Number of electrons                       :", num_electrons(1),"|"

    endif
    write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'

  end subroutine elec_report_parameters


  !=========================================================================
  subroutine elec_read_band_gradient
    !=========================================================================
    ! Read the .cst_ome file in paralell if appropriate. These are the 
    ! gradients of the bands at each kpoint.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_gradient,nspins,nbands                                       
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms, only : on_root, comms_slice, my_node_id, num_nodes, root_id,&
         & comms_recv, comms_send
    use od_io,    only : io_time, filename_len, seedname, stdout, io_file_unit,&
         & io_error  
    use od_cell,  only : num_kpoints_on_node,nkpoints
    use od_constants,only : bohr2ang, H2eV
    use od_parameters, only : legacy_file_format
    implicit none

    integer :: gradient_unit,i,ib,jb,is,ik,inodes,ierr 
    character(filename_len) :: gradient_filename

    real(kind=dp) :: time0, time1

    ! Check that we haven't already done this.

    if(allocated(band_gradient)) return 

    time0=io_time()
    if(on_root) then
       gradient_unit=io_file_unit()
       gradient_filename=trim(seedname)//".cst_ome"
       open(unit=gradient_unit,file=gradient_filename,status="old",form='unformatted',err=101)
    endif

    ! Figure out how many kpoint should be on each node
    call comms_slice(nkpoints,num_kpoints_on_node)
    allocate(band_gradient(1:nbands,1:nbands,1:3,1:num_kpoints_on_node(my_node_id),1:nspins),stat=ierr)
    if (ierr/=0) call io_error('Error: Problem allocating band_gradient in read_band_energy')

    if(legacy_file_format) then

       if(on_root) then
          do inodes=1,num_nodes-1
             do ik=1,num_kpoints_on_node(inodes)
                do is=1,nspins
                   do i=1,3
                      do jb=1,nbands
                         do ib=1,nbands
                            ! Read in units of Ha Bohr^2 / Ang
                            read (gradient_unit) band_gradient(ib,jb,i,ik,is)
                         end do
                      end do
                   end do
                end do
             end do
             call comms_send(band_gradient(1,1,1,1,1),nbands*nbands*3*nspins*num_kpoints_on_node(inodes),inodes)
          end do

          do ik=1,num_kpoints_on_node(0)
             do is=1,nspins
                do i=1,3
                   do jb=1,nbands
                      do ib=1,nbands
                         ! Read in units of Ha Bohr^2 / Ang
                         read (gradient_unit) band_gradient(ib,jb,i,ik,is)
                      end do
                   end do
                end do
             end do
          end do
       endif

    else ! sane file format

       if(on_root) then
          do inodes=1,num_nodes-1
             do ik=1,num_kpoints_on_node(inodes)
                do is=1,nspins
                   read(gradient_unit) (((band_gradient(ib,jb,i,ik,is),ib=1,nbands)&
                        ,jb=1,nbands),i=1,3)
                end do
             end do
             call comms_send(band_gradient(1,1,1,1,1),nbands*nbands*3*nspins*num_kpoints_on_node(inodes),inodes)
          end do
          do ik=1,num_kpoints_on_node(0)
             do is=1,nspins
                read(gradient_unit) (((band_gradient(ib,jb,i,ik,is),ib=1,nbands)&
                     ,jb=1,nbands),i=1,3)
             end do
          end do
       end if
    end if

    if(.not. on_root) then
       call comms_recv(band_gradient(1,1,1,1,1),nbands*nbands*3*nspins*num_kpoints_on_node(my_node_id),root_id)
    end if

    if(on_root) close (unit=gradient_unit)


    ! Convert all band gradients to eV Ang
    if(legacy_file_format) then
       band_gradient=band_gradient*bohr2ang*bohr2ang*H2eV
    else
       band_gradient=band_gradient*bohr2ang*H2eV
    end if

    time1=io_time()
    if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to read band gradients ',time1-time0,' (sec)'

    return

101 call io_error('Error: Problem opening cst_ome file in read_band_gradient')

  end subroutine elec_read_band_gradient


  !=========================================================================
  subroutine elec_read_band_energy !(band_energy,kpoint_r,kpoint_weight)
    !=========================================================================
    ! Read the .bands file in the kpoint list, kpoint weights and band energies
    ! also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_energy, efermi_castep, num_electrons      
    ! spin_polarised, electrons_per_state, nspins,nbands                             
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_constants, only : H2eV
    use od_cell,      only : nkpoints,kpoint_r,kpoint_weight,cell_find_MP_grid,&
         & real_lattice,kpoint_grid_dim,num_kpoints_on_node
    use od_comms,     only : comms_bcast,comms_send,comms_recv,num_nodes,my_node_id,&
         & on_root,root_id,comms_slice
    use od_io,        only : io_file_unit, seedname, filename_len,stdout, io_time,&
         & io_error

    implicit none

    real(kind=dp), allocatable :: all_kpoints(:,:)
    integer :: inodes,ik,is,ib,band_unit,iall_kpoints,i
    integer :: all_kpoints_pointer, dum_i1, ierr, str_pos
    character(filename_len) :: band_filename
    character(len=80) :: dummy
    real(kind=dp) :: time0, time1

    time0=io_time()

    ! Check that we haven't already read in the energies
    if(allocated(band_energy)) return

    !Open the bands file
    band_unit=io_file_unit()
    band_filename=trim(seedname)//".bands"

    ! Read the header from the bands file
    if(on_root) then
       open(unit=band_unit,file=band_filename,status="old",form='formatted',err=100)
       read (band_unit,'(a)') dummy
       str_pos=index(dummy,'k-points')
       read (dummy(str_pos+8:),*) nkpoints
       read (band_unit,'(a)') dummy
       str_pos=index(dummy,'components')
       read (dummy(str_pos+10:),*) nspins
       read (band_unit,'(a)') dummy

       allocate(num_electrons(nspins),stat=ierr)
       if(ierr/=0) stop " Error : cannot allocate num_electrons"
       str_pos=index(dummy,'electrons')
       read (dummy(str_pos+10:),*) num_electrons(:)
       read (band_unit,'(a)') dummy
       str_pos=index(dummy,'eigenvalues')
       read (dummy(str_pos+11:),*) nbands
       read (band_unit,'(a)') dummy
       str_pos=index(dummy,'units)')
       read (dummy(str_pos+6:),'(f12.4)') efermi_castep
       read (band_unit,'(a)') dummy
       read (band_unit,*) real_lattice(:,1)
       read (band_unit,*) real_lattice(:,2)
       read (band_unit,*) real_lattice(:,3)
    end if

    call comms_bcast(nspins,1)
    if(.not.on_root) then
       allocate(num_electrons(nspins),stat=ierr)
       if(ierr/=0) stop " Error : cannot allocate num_electrons"
    endif
    call comms_bcast(num_electrons(1),nspins)
    call comms_bcast(nkpoints,1)
    call comms_bcast(nbands,1)
    call comms_bcast(efermi_castep,1)
    !
    call comms_slice(nkpoints,num_kpoints_on_node)
    !
    allocate(band_energy(1:nbands,1:nspins,1:num_kpoints_on_node(my_node_id)),stat=ierr)
    if (ierr/=0) call io_error('Error: Problem allocating band_energy in read_band_energy') 
    allocate(kpoint_weight(1:num_kpoints_on_node(my_node_id)),stat=ierr)
    if (ierr/=0)  call io_error('Error: Problem allocating kpoint_weight in read_band_energy')  
    allocate(kpoint_r(1:3,1:num_kpoints_on_node(my_node_id)),stat=ierr)
    if (ierr/=0)  call io_error('Error: Problem allocating kpoint_r in read_band_energy') 

    if(on_root) then
       allocate(all_kpoints(1:3,nkpoints),stat=ierr)
       if (ierr/=0)  call io_error('Error: Problem allocating all_kpoints in read_band_energy')
       iall_kpoints=1
       do inodes=1,num_nodes-1
          do ik=1,num_kpoints_on_node(inodes)
             read (band_unit,'(a)') dummy
             str_pos=index(dummy,'K-point')
             read (dummy(str_pos+7:),*) dum_i1, kpoint_r(1,ik), kpoint_r(2,ik), kpoint_r(3,ik), kpoint_weight(ik)
             do i=1,3
                all_kpoints(i,iall_kpoints)=kpoint_r(i,ik)
             end do
             iall_kpoints=iall_kpoints+1
             do is=1,nspins
                read (band_unit,*) dummy
                do ib=1,nbands
                   read (band_unit,*) band_energy(ib,is,ik) !NB spin <-> kpt swapped
                end do
             end do
          end do
          call comms_send(band_energy(1,1,1),nbands*nspins*num_kpoints_on_node(inodes),inodes)
          call comms_send(kpoint_r(1,1),3*num_kpoints_on_node(inodes),inodes)
          call comms_send(kpoint_weight(1),num_kpoints_on_node(inodes),inodes)
       end do

       do ik=1,num_kpoints_on_node(0)
          read (band_unit,'(a)') dummy
          str_pos=index(dummy,'K-point')
          read (dummy(str_pos+7:),*) dum_i1,kpoint_r(1,ik), kpoint_r(2,ik), kpoint_r(3,ik), kpoint_weight(ik)
          do i=1,3
             all_kpoints(i,iall_kpoints)=kpoint_r(i,ik)
          enddo
          iall_kpoints=iall_kpoints+1
          do is=1,nspins
             read (band_unit,*) dummy
             do ib=1,nbands
                read (band_unit,*) band_energy(ib,is,ik) !NB spin <-> kpt swapped
             end do
          end do
       end do

       ! Do this here so we can free up the all_kpoints memory
       call cell_find_MP_grid(all_kpoints,nkpoints,kpoint_grid_dim)
       deallocate(all_kpoints,stat=ierr) 
       if (ierr/=0)  call io_error('Error: Problem deallocating all_kpoints in read_band_energy')
    endif

    if(.not. on_root) then
       call comms_recv(band_energy(1,1,1),nbands*nspins*num_kpoints_on_node(my_node_id),root_id)
       call comms_recv(kpoint_r(1,1),3*num_kpoints_on_node(my_node_id),root_id)
       call comms_recv(kpoint_weight(1),num_kpoints_on_node(my_node_id),root_id)
    end if


    if(on_root) close (unit=band_unit)

    band_energy=band_energy*H2eV
    efermi_castep=efermi_castep*H2eV

    ! Things that follow
    if(nspins .lt. 2) then
       spin_polarised=.false.
       electrons_per_state=2.0_dp
    else
       spin_polarised=.true.
       electrons_per_state=1.0_dp
    endif

    time1=io_time()
    if(on_root) write(stdout,'(1x,a40,f11.3,a)') 'Time to read band energies  ',time1-time0,' (sec)'

    return
100 call io_error('Error: Problem opening bands file in read_band_energy') 
  end subroutine elec_read_band_energy


  !=========================================================================
  subroutine elec_read_elnes_mat 
    !=========================================================================
    ! Read the .bands file in the kpoint list, kpoint weights and band energies
    ! also obtain, nkpoints, nspins, num_electrons(:),nbands, efermi_castep
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: band_energy, efermi_castep, num_electrons      
    ! spin_polarised, electrons_per_state, nspins,nbands                             
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_constants, only : H2eV
    use od_cell,      only : nkpoints,kpoint_r,kpoint_weight,cell_find_MP_grid,&
         & real_lattice,kpoint_grid_dim,num_kpoints_on_node
    use od_comms,     only : comms_bcast,comms_send,comms_recv,num_nodes,my_node_id,&
         & on_root,root_id,comms_slice
    use od_io,        only : io_file_unit, seedname, filename_len,stdout, io_time,&
         & io_error
    use od_parameters, only : legacy_file_format

    implicit none

    real(kind=dp), allocatable :: all_kpoints(:,:)
    integer :: inodes,ik,is,ib,nk,ns,nb,indx
    integer :: dum_i1, ierr, str_pos,elnes_unit,orb
    character(filename_len) :: elnes_filename
    character(len=80) :: dummy
    real(kind=dp) :: time0, time1

    time0=io_time()

    ! Check that we haven't already read in the energies
    if(allocated(elnes_mat)) return

    !Open the band sfile
    elnes_unit=io_file_unit()
    elnes_filename=trim(seedname)//".eels_mat"

    ! Read the details of the core orbitals from the elnes file
    if(on_root) then
       open(unit=elnes_unit,file=elnes_filename,form='unformatted',err=100,status='old')

       read(elnes_unit) elnes_mwab%norbitals
       read(elnes_unit) elnes_mwab%nbands   
       read(elnes_unit) elnes_mwab%nkpoints 
       read(elnes_unit) elnes_mwab%nspins   

       ! check these agree with band data?

       allocate(elnes_orbital(elnes_mwab%norbitals),stat=ierr)
       if(ierr/=0) call io_error(' Error : cannot allocate elnes_orbital')

       read(elnes_unit) elnes_orbital(1:elnes_mwab%norbitals)%species_no      
       read(elnes_unit) elnes_orbital(1:elnes_mwab%norbitals)%rank_in_species 
       read(elnes_unit) elnes_orbital(1:elnes_mwab%norbitals)%shell           
       read(elnes_unit) elnes_orbital(1:elnes_mwab%norbitals)%am_channel      

    end if
    call comms_bcast(elnes_mwab%norbitals,1)
    call comms_bcast(elnes_mwab%nbands,1)
    call comms_bcast(elnes_mwab%nkpoints,1)
    call comms_bcast(elnes_mwab%nspins,1)
    if(.not. on_root) then
       allocate(elnes_orbital(elnes_mwab%norbitals),stat=ierr)
       if(ierr/=0) call io_error(" Error : cannot allocate elnes_orbital")
    end if
    call comms_bcast(elnes_orbital(1)%species_no      ,elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital(1)%rank_in_species ,elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital(1)%shell           ,elnes_mwab%norbitals)
    call comms_bcast(elnes_orbital(1)%am_channel      ,elnes_mwab%norbitals)

    ! assume same data distribution as bands

    allocate(elnes_mat(1:elnes_mwab%norbitals,1:elnes_mwab%nbands,1:3, &
         1:num_kpoints_on_node(my_node_id),1:elnes_mwab%nspins),stat=ierr)
    if (ierr/=0) call io_error('Error: Problem allocating elnes_mat in read_band_energy')

    if(on_root) then
       if(legacy_file_format) then
          do inodes=1,num_nodes-1
             do ik=1,num_kpoints_on_node(inodes)
                do ns=1,elnes_mwab%nspins
                   do orb=1,elnes_mwab%norbitals
                      do nb=1,elnes_mwab%nbands
                         read(elnes_unit) (elnes_mat(orb,nb,indx,ik,ns),indx=1,3)
                      end do
                   end do
                end do
             end do
             call comms_send(elnes_mat(1,1,1,1,1),elnes_mwab%norbitals*elnes_mwab%nbands*3*&
                  nspins*num_kpoints_on_node(inodes),inodes)
          end do
          
          do ik=1,num_kpoints_on_node(0)
             do ns=1,elnes_mwab%nspins
                do orb=1,elnes_mwab%norbitals
                   do nb=1,elnes_mwab%nbands
                      read(elnes_unit) (elnes_mat(orb,nb,indx,ik,ns),indx=1,3)
                   end do
                end do
             end do
          end do

       else ! sane format
          do inodes=1,num_nodes-1
             do ik=1,num_kpoints_on_node(inodes)
                do ns=1,elnes_mwab%nspins
                   read(elnes_unit) (((elnes_mat(orb,nb,indx,ik,ns),orb=1,elnes_mwab%norbitals),&
                        nb=1,elnes_mwab%nbands),indx=1,3)
                end do
             end do
             call comms_send(elnes_mat(1,1,1,1,1),elnes_mwab%norbitals*elnes_mwab%nbands*3*nspins*&
                  num_kpoints_on_node(inodes),inodes)
          end do
          
          do ik=1,num_kpoints_on_node(0)
             do ns=1,elnes_mwab%nspins
                read(elnes_unit) (((elnes_mat(orb,nb,indx,ik,ns),orb=1,elnes_mwab%norbitals),&
                     nb=1,elnes_mwab%nbands),indx=1,3)
             end do
          end do
       end if
    endif

    if(.not. on_root) then
       call comms_recv(elnes_mat(1,1,1,1,1),elnes_mwab%norbitals*elnes_mwab%nbands*3*nspins*&
            num_kpoints_on_node(my_node_id),root_id)
    end if
    
    if(on_root) close(elnes_unit)
    
    return

100 call io_error('Error: Problem opening elnes file in elec_read_elnes_mat') 


  end subroutine elec_read_elnes_mat

  !=========================================================================
  subroutine elec_pdos_read
    !=========================================================================
    ! Read in the full pdos_weights. Write out any variables that we find on the 
    ! way. These will be checked for consistency in the dos module. We can't do it 
    ! yet as we haven't read the bands file.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: pw, pdos_weights, pdos_orbital                          
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms,     only : on_root
    use od_io,        only : io_file_unit, io_error, seedname, stdout

    implicit none

    ! Band indices used in the read-in of the pdos 
    integer, allocatable, dimension(:,:) :: nbands_occ
    real(kind=dp)                        :: dummyr1,dummyr2,dummyr3
    integer                              :: dummyi,ib,ik,is
    integer                              :: pdos_in_unit,ios,ierr

    pdos_in_unit=io_file_unit()

    if(allocated(pdos_orbital)) deallocate(pdos_orbital)

    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R
    open (pdos_in_unit, iostat=ios, status='old', file=trim(seedname)//".pdos_weights", form='unformatted')
    if(ios.ne.0) call io_error ("Error : Cannot open pDOS weights")

    read(pdos_in_unit) pdos_mwab%nkpoints
    read(pdos_in_unit) pdos_mwab%nspins
    read(pdos_in_unit) pdos_mwab%norbitals
    read(pdos_in_unit) pdos_mwab%nbands

    allocate(pdos_orbital(pdos_mwab%norbitals),stat=ierr)
    if(ierr/=0) stop " Error : cannot allocate orbital"

    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%species_no
    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%rank_in_species
    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%am_channel
    !-------------------------------------------------------------------------!


    !-------------------------------------------------------------------------!
    ! N O W   R E A D   T H E   D A T A
    allocate(nbands_occ(1:pdos_mwab%nkpoints,1:pdos_mwab%nspins),stat=ierr)
    if(ierr/=0) stop " Error : cannot allocate nbands_occ"
    allocate(pdos_weights(1:pdos_mwab%norbitals,1:pdos_mwab%nbands,1:pdos_mwab%nkpoints,1:pdos_mwab%nspins),stat=ierr)
    if(ierr/=0) stop " Error : cannot allocate pdos_weights"

    do ik=1,pdos_mwab%nkpoints
       ! The kpoint number, followed by the kpoint-vector
       read(pdos_in_unit) dummyi, dummyr1, dummyr2, dummyr3
       do is=1, pdos_mwab%nspins
          read(pdos_in_unit) dummyi ! this is the spin number
          read(pdos_in_unit) nbands_occ(ik,is)
          do ib=1,nbands_occ(ik,is)
             read(pdos_in_unit) pdos_weights(1:pdos_mwab%norbitals,ib,ik,is)
          enddo
       enddo
    enddo
    !-------------------------------------------------------------------------!


    !-------------------------------------------------------------------------!
    ! O U T P U T   O U R   F I N D I N G S
    ! call count_atoms(pdos_orbital,pdos_mwab%norbitals,nions)

    if(on_root) then
       write(stdout, *)
       write(stdout,'(1x,a78)') '+------------------------- Partial DOS Band Data ----------------------------+'
       write(stdout,'(1x,a46,i4,a28)')   '|  Number of Bands                           :',pdos_mwab%nbands,"|"
       write(stdout,'(1x,a46,i4,a28)')   '|  Number of K-points                        :',pdos_mwab%nkpoints,"|"
       write(stdout,'(1x,a46,i4,a28)')   '|  Number of LCAO                            :',pdos_mwab%norbitals,"|"
       if(pdos_mwab%nspins > 1) then
          write(stdout,'(1x,a78)') '|  Spin-Polarised Calculation                :  True                         |'
       else
          write(stdout,'(1x,a78)') '|  Spin-Polarised Calculation                :  False                        |'
       endif
       !     write(stdout,'(1x,a46,i4,a28)')   '|  Number of Ions                            :',nions,"|"
       write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
    endif
    !-------------------------------------------------------------------------!


    !-------------------------------------------------------------------------!
    ! F I N A L I S E   
    close(unit=pdos_in_unit)
    if (allocated(nbands_occ)) deallocate(nbands_occ,stat=ierr)
    if (ierr/=0) stop " Error: cannot deallocate nbands_occ"
    !-------------------------------------------------------------------------!
  end subroutine elec_pdos_read

  !=========================================================================
  subroutine elec_pdos_read_orbitals
    !=========================================================================
    ! Read in the full pdos_weights. Write out any variables that we find on the 
    ! way. These will be checked for consistency in the dos module. We can't do it 
    ! yet as we haven't read the bands file.
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables: pw, pdos_weights, pdos_orbital                          
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by  A J Morris                                         Dec 2010
    !=========================================================================
    use od_comms,     only : on_root
    use od_io,        only : io_file_unit, io_error, seedname, stdout

    implicit none

    ! Band indices used in the read-in of the pdos 
    integer, allocatable, dimension(:,:) :: nbands_occ
    real(kind=dp)                        :: dummyr1,dummyr2,dummyr3
    integer                              :: dummyi,ib,ik,is
    integer                              :: pdos_in_unit,ios,ierr

    pdos_in_unit=io_file_unit()

    !-------------------------------------------------------------------------!
    ! R E A D   T H E   D A T A   H E A D E R
    open (pdos_in_unit, iostat=ios, status='old', file=trim(seedname)//".pdos_weights", form='unformatted')
    if(ios.ne.0) call io_error ("Error : Cannot open pDOS weights")

    read(pdos_in_unit) pdos_mwab%nkpoints
    read(pdos_in_unit) pdos_mwab%nspins
    read(pdos_in_unit) pdos_mwab%norbitals
    read(pdos_in_unit) pdos_mwab%nbands

    allocate(pdos_orbital(pdos_mwab%norbitals),stat=ierr)
    if(ierr/=0) stop " Error : cannot allocate orbital"

    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%species_no
    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%rank_in_species
    read(pdos_in_unit) pdos_orbital(1:pdos_mwab%norbitals)%am_channel
    !-------------------------------------------------------------------------!


    ! F I N A L I S E   
    close(unit=pdos_in_unit)

    !-------------------------------------------------------------------------!
  end subroutine elec_pdos_read_orbitals


  subroutine elec_dealloc_pdos
    use od_io, only : io_error
    implicit none
    integer :: ierr

    if(allocated(pdos_weights)) then
       deallocate(pdos_weights,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating pdos_weights in elec_dealloc_pdos')
    end if

    if(allocated(pdos_orbital)) then
       deallocate(pdos_orbital,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating pdos_orbital in elec_dealloc_pdos')
    end if



  end subroutine elec_dealloc_pdos

  subroutine elec_dealloc_elnes
    use od_io, only : io_error
    implicit none
    integer :: ierr

    if(allocated(elnes_mat)) then
       deallocate(elnes_mat,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating elnes_mat in elec_dealloc_elnes')
    end if

    if(allocated(elnes_orbital)) then
       deallocate(elnes_orbital,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating elnes_orbital in elec_dealloc_elnes')
    end if


  end subroutine elec_dealloc_elnes

  subroutine elec_dealloc_band_gradient
    use od_io, only : io_error
    implicit none
    integer :: ierr

    if(allocated(band_gradient)) then
       deallocate(band_gradient,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating band_gradient in elec_dealloc_band_gradient')
    end if

  end subroutine elec_dealloc_band_gradient



end module od_electronic
