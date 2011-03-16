!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
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
  ! use od_dos,       only : matrix_weights_array_boundaries
  ! use od_cell,      only : nkpoints
  ! use od_comms,     only : on_root
  ! use od_io,        only : io_file_unit, io_error, seedname, stdout
  ! use od_dos,       only : E, matrix_weights_array_boundaries
  ! use od_electronic,only : nbands,nspins
  ! use od_parameters

  implicit none


  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  real(kind=dp), public, allocatable, save  :: dos_partial(:,:,:)
  !-------------------------------------------------------------------------!

  private

  ! public :: pdos_read
  public :: pdos_merge
  public :: pdos_write

  ! integer :: nions ! This will probably go into the ion module.

contains



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
    use od_dos_utils,       only : E
    use od_parameters,only : nbins
    use od_electronic, only         : pdos_mwab

    implicit none

    character(len=20) :: string 
    integer :: idos, i
    write(string,'(I4,"(x,es14.7)")') pdos_mwab%norbitals

    if(pdos_mwab%nspins>1) then
       dos_partial(:,2,:)=-dos_partial(:,2,:)
       do idos = 1,nbins
          write(58,'(es14.7,'//trim(string)//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=1,pdos_mwab%norbitals) &
               & ,(dos_partial(idos,2,i),i=1,pdos_mwab%norbitals)
       end do
    else
       do idos = 1,nbins
          write(58,'(es14.7,'//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=1,pdos_mwab%norbitals)
       end do
    endif


  end subroutine pdos_write
!!$
!!$!===============================================================================
!!$ subroutine count_atoms(orbital,num_orbitals,num_atoms)
!!$!===============================================================================
!!$! From the program LinDOS (AJM)
!!$! Take the orbial information and work out the number of atoms that the LCAO 
!!$! describe
!!$!===============================================================================
!!$   use od_io, only : io_error
!!$   implicit none
!!$   integer, intent(in)            :: num_orbitals
!!$   type(orbitals), intent(inout)  :: orbital(1:num_orbitals) ! sepcies, ! num of each species ! l channel
!!$   integer, intent(out)           :: num_atoms
!!$
!!$   integer, allocatable           :: species_count(:)
!!$   integer                        :: num_species, ion_count
!!$   integer                        :: i, ierr
!!$
!!$   num_species=maxval(orbital(:)%species_no)  ! The maximum value is the highest species rank
!!$
!!$   allocate(species_count(1:num_species), stat=ierr)
!!$   if(ierr/=0) call io_error( " Error : cannot allocate species_count")
!!$
!!$   species_count=0
!!$   ion_count=0
!!$
!!$   do i=1,num_orbitals
!!$    ! If the species number is greater than the number we have for that species then use this 
!!$    ! new number instead
!!$    ! NB I'm using data from the array orbital to index species count! :S
!!$    if(orbital(i)%rank_in_species>species_count(orbital(i)%species_no)) then
!!$      species_count(orbital(i)%species_no)=orbital(i)%rank_in_species
!!$      ion_count=ion_count+1
!!$    endif
!!$    orbital(i)%ion_no=ion_count
!!$   enddo
!!$
!!$   num_atoms=sum(species_count(:))
!!$   
!!$   if(allocated(species_count)) then
!!$     deallocate(species_count, stat=ierr)   
!!$     if(ierr/=0) stop " Error : cannot deallocate  species_count"
!!$   endif
!!$  end subroutine count_atoms

end module od_pdos
