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
  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights

  !-------------------------------------------------------------------------!

  private

  public :: pdos_calculate

  ! integer :: nions ! This will probably go into the ion module.

contains


  subroutine pdos_calculate
    use od_electronic, only :  pdos_orbital, pdos_weights,pdos_mwab
    use od_dos_utils, only : dos_utils_calculate
    use od_comms, only : on_root

    implicit none

    call pdos_merge

    call dos_utils_calculate(matrix_weights, dos_partial)

    if (on_root) then
       call pdos_write
    endif



  end subroutine pdos_calculate


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
  subroutine pdos_get_string
    !===============================================================================
    ! This is a mindbendingly horrific exercise in book-keeping
    !===============================================================================
    use od_parameters, only : pdos_string
    use od_io, only : maxlen, io_error
    implicit none
    character(len=maxlen) :: ctemp, ctemp2,c_am,m_string

    integer :: pos_l,pos_r,ia,iz,idiff,ic1,ic2,num_proj,max_site,max_atoms,species,num_sites,num_am
    character(len=2)  :: c_symbol
    logical :: pdos_sum,am_sum,site_sum

    integer   :: kl, in,loop,num1,num2,i_punc,pos3,loop_l,loop_a,loop_p,max_am,num_species
    integer   :: counter,i_digit,loop_r,range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit="0123456789"
    character(len=2) , parameter :: c_range="-:"
    character(len=3) , parameter :: c_sep=" ,"
    character(len=5) , parameter :: c_punc=" ,-:"
    character(len=5)  :: c_num1,c_num2
    integer :: pdos_atoms(20),pdos_ang(4)
    integer, allocatable :: pdos_array(:,:,:,:)
    

    pdos_atoms=0;pdos_ang=0
    site_sum=.false.
    am_sum=.false.

    ctemp=pdos_string
    pdos_sum=.false.
    if(index(ctemp,'sum:')==1) then
       pdos_sum=.true.
       ctemp=ctemp(5:)
    end if
    write(*,*) 'pdos_sum',pdos_sum
    ! take 1st part of string

    ! look for Ang Mtm string eg (s,p)
    c_am=''
    pos_l=index(ctemp,'(')
    if(pos_l>0) then
       pos_r=index(ctemp,')')
       if(pos_r==0) call io_error ('found ( but no )')
       if(pos_r<=pos_l) call io_error (' ) before (')
       c_am=ctemp(pos_l+1:pos_r-1)
       write(*,*) c_am
       ctemp=ctemp(:pos_l-1)
       write(*,*) ctemp
    else ! implicit sum over AM
       am_sum=.true.
    end if

    ia = ichar('a')
    iz = ichar('z')
    idiff = ichar('Z')-ichar('z')
    
    ic1=ichar(ctemp(1:1))
    if(ic1<ia .or. ic1>iz) call io_error ('problem reading atomic symbol in pdos string')
    ic2=ichar(ctemp(2:2))
    if(ic2>=ia .and. ic1<=iz) then
       c_symbol=ctemp(1:2)
       ctemp=ctemp(3:)
    else
       c_symbol=ctemp(1:1)
       ctemp=ctemp(2:)
    end if
    write(*,*) 'c_symbol: ',c_symbol
    write(*,*) ctemp


    !Count atoms numbers
    counter=0
    dummy=adjustl(ctemp)
    if (len_trim(dummy)>0) then
       dummy=adjustl(dummy)
       do 
          i_punc=scan(dummy,c_punc)
          if(i_punc==0) call io_error('Error parsing keyword ') 
          c_num1=dummy(1:i_punc-1)
          read(c_num1,*,err=101,end=101) num1
          dummy=adjustl(dummy(i_punc:))
          !look for range
          if(scan(dummy,c_range)==1) then
             i_digit=scan(dummy,c_digit)
             dummy=adjustl(dummy(i_digit:))
             i_punc=scan(dummy,c_punc)
             c_num2=dummy(1:i_punc-1)
             read(c_num2,*,err=101,end=101) num2
             dummy=adjustl(dummy(i_punc:))
             range_size=abs(num2-num1)+1
             do loop_r=1,range_size
                counter=counter+1
                pdos_atoms(min(num1,num2)+loop_r-1)=1
             end do
          else
             counter=counter+1 
             pdos_atoms(num1)=1
          end if
          
          if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
          if(scan(dummy,c_range)==1) call io_error('Error parsing keyword incorrect range') 
          if(index(dummy,' ')==1) exit
       end do
    else
       site_sum=.true.
    end if
    
    do loop=1,20
       write(*,*) loop, pdos_atoms(loop)
    end do

    ! count am
    counter=0
    dummy=adjustl(c_am)
    print*,dummy
    if (len_trim(dummy)>0) then
       do
          pos3=index(dummy,',')
          if (pos3==0) then
             ctemp2=dummy
          else
             ctemp2=dummy(:pos3-1)
          endif
          read(ctemp2(1:),*,err=106,end=106) m_string
          select case (trim(adjustl(m_string)))
          case ('s')
             pdos_ang(1)=1
          case ('p')
             pdos_ang(2)=1
          case ('d')
             pdos_ang(3)=1
          case ('f')
             pdos_ang(4)=1
          case default
             call io_error('param_get_projection: Problem reading l state ')
          end select
          if (pos3==0) exit
          dummy=dummy(pos3+1:)
       enddo
    else
       am_sum=.true.
    endif
    
    do loop=1,4
       write(*,*) loop,pdos_ang(loop)
    end do

    if(site_sum) then 
       num_sites=1
    else
       num_sites=count(pdos_atoms==1)
    end if
    if(am_sum) then
       num_am=1
    else
       num_am=count(pdos_ang==1)
    end if
    write(*,*) 'site_sum',site_sum
    write(*,*) 'am_sum',am_sum
    num_proj=num_am*num_sites
    write(*,*) 'num_proj',num_proj

    num_species=1;    max_atoms=20;max_am=4
    allocate(pdos_array(num_species,max_atoms,max_am,num_proj))
    pdos_array=0

    species=1
    loop_p=1
    if(site_sum.and.am_sum) then
       pdos_array(species,:,:,loop_p)=1
    elseif(site_sum.and..not.am_sum) then
       do loop_l=1,max_am
          if(pdos_ang(loop_l)==0) cycle 
          pdos_array(species,:,loop_l,loop_p)=1
          loop_p=loop_p+1
       end do
    elseif(.not.site_sum.and.am_sum) then
       do loop_a=1,max_atoms
          write(*,*) loop_a
          if(pdos_atoms(loop_a)==0) cycle 
          pdos_array(species,loop_a,:,loop_p)=1
          loop_p=loop_p+1
       end do
    else
       do loop_l=1,max_am
          if(pdos_ang(loop_l)==0) cycle 
          do loop_a=1,max_atoms
             if(pdos_atoms(loop_a)==0) cycle 
             pdos_array(species,loop_a,loop_l,loop_p)=1
             loop_p=loop_p+1
          end do
       end do
    end if

    do loop=1,num_proj
       write(*,*) 'projection',loop
       write(*,*) pdos_array(:,:,:,loop)
    end do
    write(*,*) 'finish'
    


    return

101 call io_error('Error parsing keyword ')
106 call io_error('param_get_projection: Problem reading m state into string ')


  end subroutine pdos_get_string





  !===============================================================================
  subroutine pdos_write
    !===============================================================================
    ! Write out the pdos that was requested. Make a pretty header so that the user
    ! knows what each column means
    !===============================================================================
    use od_dos_utils,       only : E
    use od_parameters,only : dos_nbins
    use od_electronic, only         : pdos_mwab

    implicit none

    character(len=20) :: string 
    integer :: idos, i
    write(string,'(I4,"(x,es14.7)")') pdos_mwab%norbitals

    if(pdos_mwab%nspins>1) then
       dos_partial(:,2,:)=-dos_partial(:,2,:)
       do idos = 1,dos_nbins
          write(58,'(es14.7,'//trim(string)//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=1,pdos_mwab%norbitals) &
               & ,(dos_partial(idos,2,i),i=1,pdos_mwab%norbitals)
       end do
    else
       do idos = 1,dos_nbins
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
