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
 real(kind=dp), allocatable, public, save  :: band_gradient(:,:,:,:,:)

 real(kind=dp), public, save :: efermi ! The fermi energy we finally decide on
 real(kind=dp), public, save :: efermi_castep ! Fermi energy as reported by CASTEP

 real(kind=dp), allocatable, public, save :: num_electrons(:) ! Holds up-spin and 
                                                              ! down-spin

 integer, public, save       :: nbands, nspins
 real(kind=dp), public, save :: electrons_per_state ! 2 for non-spin-P
                                                    ! 1 for spin-P
 logical, public, save       :: spin_polarised
  !-------------------------------------------------------------------------!

 private
 
 !-------------------------------------------------------------------------!
 ! G L O B A L   F U N C T I O N S
 public :: elec_report_parameters
 public :: elec_read_band_energy
 public :: elec_read_band_gradient
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
  implicit none

  integer :: gradient_unit,i,ib,jb,is,ik,inodes,ierr 
  character(filename_len) :: gradient_filename

  real(kind=dp) :: time0, time1

  ! Check that we haven't already done this.
  if(allocated(band_gradient)) return 

  time0=io_time()
  gradient_unit=io_file_unit()
  gradient_filename=trim(seedname)//".cst_ome"

  if(on_root) then
    open(unit=gradient_unit,file=gradient_filename,status="old",form='unformatted',err=101)
  endif

 ! Figure out how many kpoint should be on each node
  call comms_slice(nkpoints,num_kpoints_on_node)

  allocate(band_gradient(1:nbands,1:nbands,1:3,1:num_kpoints_on_node(my_node_id),1:nspins),stat=ierr)
  if (ierr/=0) call io_error('Error: Problem allocating band_gradient in read_band_energy')

  if(on_root) then
    do inodes=1,num_nodes-1
      do ik=1,num_kpoints_on_node(inodes)
        do is=1,nspins
          do i=1,3
            do jb=1,nbands
              do ib=1,nbands
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
                read (gradient_unit) band_gradient(ib,jb,i,ik,is)
              end do
           end do
         end do
       end do 
   end do
  endif

  if(.not. on_root) then
     call comms_recv(band_gradient(1,1,1,1,1),nbands*nbands*3*nspins*num_kpoints_on_node(my_node_id),root_id)
  end if

  close (unit=gradient_unit)

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
  
!Open the band sfile
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

    call comms_bcast(nspins,1)
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

  call comms_bcast(nkpoints,1)
  call comms_bcast(nbands,1)
  call comms_bcast(efermi_castep,1)

  ! Broadcast the number of kpoints
  call comms_bcast(nkpoints,1)

  call comms_slice(nkpoints,num_kpoints_on_node)

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

  ! Do this here so we can free up the all_kpoints memeory
   call cell_find_MP_grid(all_kpoints,nkpoints,kpoint_grid_dim)
   deallocate(all_kpoints,stat=ierr) 
   if (ierr/=0)  call io_error('Error: Problem deallocating all_kpoints in read_band_energy')
 endif

  if(.not. on_root) then
     call comms_recv(band_energy(1,1,1),nbands*nspins*num_kpoints_on_node(my_node_id),root_id)
     call comms_recv(kpoint_r(1,1),3*num_kpoints_on_node(my_node_id),root_id)
     call comms_recv(kpoint_weight(1),num_kpoints_on_node(my_node_id),root_id)
  end if


  close (unit=band_unit)
  
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
  write(stdout,'(1x,a40,f11.3,a)') 'Time to read band energies  ',time1-time0,' (sec)'

  return
100 call io_error('Error: Problem opening bands file in read_band_energy') 
end subroutine elec_read_band_energy
 
end module od_electronic
