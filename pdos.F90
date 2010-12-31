module od_pdos
!===============================================================================
! OptaDOS is written in such a way, that the dos routines take an optinal
! matrix weights variable.  This can be a pdos weighting or any other dos weighting
! that the OptaDOS programmer may wish. 
!   To this end, pdos is just an example of a dos weighting. This pdos module 
! performs the 3 functinoals necessary. Namely,
! (1) The reading in of the pdos weights from the .pdos_weights file
! (2) The exercise in book-keeping necessary to merge and discard the weights 
!     such that the requested channels are given to the dos module
! (3) The writing out of the pdos, and labelling the columns for easy reading
!===============================================================================
 use od_constants, only : dp
 use od_dos,       only : matrix_weights_array_boundaries
! use od_cell,      only : nkpoints
! use od_comms,     only : on_root
! use od_io,        only : io_file_unit, io_error, seedname, stdout
! use od_dos,       only : E, matrix_weights_array_boundaries
! use od_electronic,only : nbands,nspins
! use od_parameters

 implicit none

!-------------------------------------------------------------------------!
! P R O T O T Y P E S
!-------------------------------------------------------------------------!
  type orbitals
   integer          :: ion_no          ! Unique ion number
   integer          :: species_no      ! Uniqie species number
   integer          :: rank_in_species ! Unique ion number within species
   integer          :: am_channel      ! The angular momentum Channel (l)
   logical          :: calc_pdos       ! Should pDoS be calculated for this 
                                       ! orbital?                                         
   character(len=1) :: am_channel_name ! Name of angular momentum channel, 
                                       ! s,p,d, etc
  end type orbitals
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! G L O B A L   V A R I A B L E S
!-------------------------------------------------------------------------!
 type(orbitals), public, allocatable, save :: pdos_orbital(:)

 real(kind=dp), public, allocatable, save  :: pdos_weights(:,:,:,:)
 real(kind=dp), public, allocatable, save  :: dos_partial(:,:,:)
!-------------------------------------------------------------------------!

 private

 public :: pdos_read
 public :: pdos_merge
 public :: pdos_write

 integer :: nions ! This will probably go into the ion module.
 
 type(matrix_weights_array_boundaries) :: pw

 contains

!=========================================================================
 subroutine pdos_read
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

  read(pdos_in_unit) pw%nkpoints
  read(pdos_in_unit) pw%nspins
  read(pdos_in_unit) pw%norbitals
  read(pdos_in_unit) pw%nbands

  allocate(pdos_orbital(pw%norbitals),stat=ierr)
  if(ierr/=0) stop " Error : cannot allocate orbital"

  read(pdos_in_unit) pdos_orbital(1:pw%norbitals)%species_no
  read(pdos_in_unit) pdos_orbital(1:pw%norbitals)%rank_in_species
  read(pdos_in_unit) pdos_orbital(1:pw%norbitals)%am_channel
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! N O W   R E A D   T H E   D A T A
 allocate(nbands_occ(1:pw%nkpoints,1:pw%nspins),stat=ierr)
 if(ierr/=0) stop " Error : cannot allocate nbands_occ"
 allocate(pdos_weights(1:pw%norbitals,1:pw%nbands,1:pw%nkpoints,1:pw%nspins),stat=ierr)
 if(ierr/=0) stop " Error : cannot allocate pdos_weights"

 do ik=1,pw%nkpoints
 ! The kpoint number, followed by the kpoint-vector
    read(pdos_in_unit) dummyi, dummyr1, dummyr2, dummyr3
   do is=1, pw%nspins
     read(pdos_in_unit) dummyi ! this is the spin number
     read(pdos_in_unit) nbands_occ(ik,is)
     do ib=1,nbands_occ(ik,is)
       read(pdos_in_unit) pdos_weights(1:pw%norbitals,ib,ik,is)
     enddo
   enddo
 enddo
!-------------------------------------------------------------------------!


!-------------------------------------------------------------------------!
! O U T P U T   O U R   F I N D I N G S
 call count_atoms(pdos_orbital,pw%norbitals,nions)

  if(on_root) then
    write(stdout, *)
    write(stdout,'(1x,a78)') '+------------------------- Partial DOS Band Data ----------------------------+'
    write(stdout,'(1x,a46,i4,a28)')   '|  Number of Bands                           :',pw%nbands,"|"
    write(stdout,'(1x,a46,i4,a28)')   '|  Number of K-points                        :',pw%nkpoints,"|"
    write(stdout,'(1x,a46,i4,a28)')   '|  Number of LCAO                            :',pw%norbitals,"|"
    if(pw%nspins > 1) then
      write(stdout,'(1x,a78)') '|  Spin-Polarised Calculation                :  True                         |'
    else
      write(stdout,'(1x,a78)') '|  Spin-Polarised Calculation                :  False                        |'
     endif
     write(stdout,'(1x,a46,i4,a28)')   '|  Number of Ions                            :',nions,"|"
     write(stdout,'(1x,a78)') '+----------------------------------------------------------------------------+'
  endif
!-------------------------------------------------------------------------!
 

!-------------------------------------------------------------------------!
! F I N A L I S E   
  close(unit=pdos_in_unit)
  if (allocated(nbands_occ)) deallocate(nbands_occ,stat=ierr)
  if (ierr/=0) stop " Error: cannot deallocate nbands_occ"
!-------------------------------------------------------------------------!
 end subroutine pdos_read

 
!===============================================================================
subroutine pdos_merge
!===============================================================================
! This is a mindbendingly horrific exercise in book-keeping
!===============================================================================
 implicit none
 ! Good isn't it?

 return
end subroutine pdos_merge


!===============================================================================
 subroutine pdos_write
!===============================================================================
! Write out the pdos that was requested. Make a pretty header so that the user
! knows what each column means
!===============================================================================
  use od_dos,       only : E
  use od_parameters,only : nbins

  implicit none

  character(len=20) :: string 
  integer :: idos, i
  write(string,'(I4,"(x,es14.7)")') pw%norbitals

  if(pw%nspins>1) then
     dos_partial(:,2,:)=-dos_partial(:,2,:)
     do idos = 1,nbins
       write(58,'(es14.7,'//trim(string)//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=1,pw%norbitals) &
          & ,(dos_partial(idos,2,i),i=1,pw%norbitals)
     end do
   else
     do idos = 1,nbins
       write(58,'(es14.7,'//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=1,pw%norbitals)
     end do
   endif

  
 end subroutine pdos_write

!===============================================================================
 subroutine count_atoms(orbital,num_orbitals,num_atoms)
!===============================================================================
! From the program LinDOS (AJM)
! Take the orbial information and work out the number of atoms that the LCAO 
! describe
!===============================================================================
   use od_io, only : io_error
   implicit none
   integer, intent(in)            :: num_orbitals
   type(orbitals), intent(inout)  :: orbital(1:num_orbitals) ! sepcies, ! num of each species ! l channel
   integer, intent(out)           :: num_atoms

   integer, allocatable           :: species_count(:)
   integer                        :: num_species, ion_count
   integer                        :: i, ierr

   num_species=maxval(orbital(:)%species_no)  ! The maximum value is the highest species rank

   allocate(species_count(1:num_species), stat=ierr)
   if(ierr/=0) call io_error( " Error : cannot allocate species_count")

   species_count=0
   ion_count=0

   do i=1,num_orbitals
    ! If the species number is greater than the number we have for that species then use this 
    ! new number instead
    ! NB I'm using data from the array orbital to index species count! :S
    if(orbital(i)%rank_in_species>species_count(orbital(i)%species_no)) then
      species_count(orbital(i)%species_no)=orbital(i)%rank_in_species
      ion_count=ion_count+1
    endif
    orbital(i)%ion_no=ion_count
   enddo

   num_atoms=sum(species_count(:))
   
   if(allocated(species_count)) then
     deallocate(species_count, stat=ierr)   
     if(ierr/=0) stop " Error : cannot deallocate  species_count"
   endif
  end subroutine count_atoms

end module od_pdos
