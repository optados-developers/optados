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

  implicit none


  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  real(kind=dp), public, allocatable, save  :: dos_partial(:,:,:)
  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights
  integer, allocatable :: pdos_projection_array(:,:,:,:)
  integer :: num_proj

  integer, parameter :: max_am=4                   ! s,p,d,f  hard coded!

  ! Data derived from the info in the pdos_weights file 
  character(len=3), allocatable :: pdos_symbol(:)   ! symbols
  integer, allocatable :: pdos_am(:,:)              ! angular mtm (num_species,max_am)
  integer, allocatable :: pdos_sites(:)             ! number of each species
  logical :: shortcut


  !-------------------------------------------------------------------------!

  private

  public :: pdos_calculate

contains


  subroutine pdos_calculate
    use od_electronic, only :elec_pdos_read
    use od_dos_utils, only : dos_utils_calculate
    use od_comms, only : on_root
    use od_parameters, only : iprint
    use od_io, only : stdout

    implicit none

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+=============== Partial Density Of States Calculation ======================+'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif

    ! read in the pdos weights
    call elec_pdos_read

    ! look at the orbitals and figure out which atoms / states we have
    call pdos_analyse_orbitals

    ! parse the pdos string to see what we want
    call pdos_get_string

    ! form the right matrix elements
    call pdos_merge

    if(on_root.and.(iprint>2)) then
       call pdos_report_projectors
    endif

    ! now compute the weighted dos
    call dos_utils_calculate(matrix_weights, dos_partial)

    ! and write everything out
    if (on_root) then
       call pdos_write
    endif

   if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+============== Partial Density Of States Calculation End ===================+'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif

  end subroutine pdos_calculate


  !===============================================================================
  subroutine pdos_merge
    !===============================================================================
    ! This is a mindbendingly horrific exercise in book-keeping
    !===============================================================================
    use od_electronic, only :  pdos_orbital, pdos_weights,pdos_mwab,nspins
    use od_cell, only : num_kpoints_on_node
    use od_comms, only : my_node_id
    use od_io, only : io_error
    implicit none

    integer :: N,N_spin,n_eigen,nproj,orb,ierr

    allocate(matrix_weights(num_proj,pdos_mwab%nbands,num_kpoints_on_node(my_node_id),nspins),stat=ierr)
    if(ierr/=0) call io_error('Error: pdos_merge - allocation of matrix_weights failed')
    matrix_weights=0.0_dp

    do N=1,num_kpoints_on_node(my_node_id) ! Loop over kpoints
       do N_spin=1,nspins                  ! Loop over spins
          do n_eigen=1,pdos_mwab%nbands    ! Loop over unoccupied states
             do nproj=1,num_proj
                do orb=1,pdos_mwab%norbitals
                   if(pdos_projection_array(pdos_orbital%species_no(orb),pdos_orbital%rank_in_species(orb) &
                        ,pdos_orbital%am_channel(orb)+1,nproj)==1) then
                      matrix_weights(nproj,n_eigen,N,N_spin)=matrix_weights(nproj,n_eigen,N,N_spin)+&
                           pdos_weights(orb,n_eigen,N,N_spin)
                   end if
                end do
             end do
          end do
       end do
    end do

    return
  end subroutine pdos_merge



  !===============================================================================
  subroutine pdos_get_string
    !===============================================================================
    ! This is a mindbendingly horrific exercise in book-keeping
    !===============================================================================
    use od_parameters, only : pdos_string
    use od_cell, only : num_species,atoms_species_num
    use od_io, only : maxlen,io_error
    implicit none
    character(len=maxlen) :: ctemp, ctemp2,ctemp3


    integer :: loop4,loop3,loop2,ierr
    logical :: pdos_sum

    integer   :: loop,pos,loop_l,loop_a,loop_p
    integer   :: i_digit,species_count,species_proj
    character(len=1) , parameter :: c_sep=":"
    integer, allocatable :: pdos_temp(:,:,:,:)

    !Check for any short cuts
    shortcut=.false.

    ctemp=pdos_string
    if(index(ctemp,'species_ang')>0) then
       num_proj=0
       do loop=1,num_species
          num_proj=num_proj+count(pdos_am(loop,:)==1)
       end do
       allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')

       pdos_projection_array=0
       loop_p=1
       do loop=1,num_species
          do loop_l=1,max_am
             if(pdos_am(loop,loop_l)==0) cycle 
             pdos_projection_array(loop,:,loop_l,loop_p)=1
             loop_p=loop_p+1
          end do
       end do
       shortcut=.true.
    elseif(index(ctemp,'species')>0) then
       num_proj=num_species
       allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')

       pdos_projection_array=0
       do loop=1,num_species
          pdos_projection_array(loop,:,:,loop)=1
       end do
       shortcut=.true.
    elseif(index(ctemp,'sites')>0) then
       num_proj=0
       do loop=1,num_species
          num_proj=num_proj+pdos_sites(loop)
       end do
       allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')

       pdos_projection_array=0
       loop_p=1
       do loop=1,num_species
          do loop_a=1,pdos_sites(loop)
             pdos_projection_array(loop,loop_a,:,loop_p)=1
             loop_p=loop_p+1
          end do
       end do
       shortcut=.true.
    elseif(index(ctemp,'angular')>0) then
       num_proj=max_am
       allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')
       pdos_projection_array=0
       loop_p=0
       do loop=1,num_proj
          pdos_projection_array(:,:,loop,loop)=1
       end do
       shortcut=.true.
    endif


    if(.not. shortcut) then

       !look for sum
       ctemp=pdos_string
       pdos_sum=.false.
       if(index(ctemp,'sum:')==1) then
          pdos_sum=.true.
          ctemp=ctemp(5:)
       end if
       ! take 1st part of string

       ctemp2=ctemp
       species_count=1;num_proj=0
       do
          !look for each species section
          ! and pass to find number of projections
          pos=index(ctemp2,c_sep)
          if(pos==0) then
             ctemp3=ctemp2
          else
             ctemp3=ctemp2(1:pos-1)
          end if
          call pdos_analyse_substring(ctemp3,species_proj)
          num_proj=num_proj+species_proj
          if (pos==0) exit
          species_count=species_count+1
          ctemp2=ctemp2(pos+1:)
       end do

       ! now allocate the correct sized array

       allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')
       pdos_projection_array=0

       ctemp2=ctemp
       do loop=1,species_count
          !loop for each species section
          !and pass fill in projection
          pos=index(ctemp2,c_sep)
          if(pos==0) then
             ctemp3=ctemp2
          else
             ctemp3=ctemp2(1:pos-1)
          end if
          call pdos_analyse_substring(ctemp3)
          ctemp2=ctemp2(pos+1:)
       end do

       if(pdos_sum) then
          allocate(pdos_temp(num_species,maxval(atoms_species_num),max_am,1),stat=ierr)
       if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_temp failed')

          pdos_temp=0
          do loop4=1,num_proj
             do loop3=1,max_am
                do loop2=1,maxval(atoms_species_num)
                   do loop=1,num_species
                      if(pdos_projection_array(loop,loop2,loop3,loop4)==1) then
                         pdos_temp(loop,loop2,loop3,1)=1
                      end if
                   end do
                end do
             end do
          end do
          deallocate(pdos_projection_array,stat=ierr)
          if(ierr/=0) call io_error('Error: pdos_get_string - deallocation of pdos_projection_array failed')
          num_proj=1
          allocate(pdos_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
          if(ierr/=0) call io_error('Error: pdos_get_string - allocation of pdos_projection_array failed')
          pdos_projection_array=0
          pdos_projection_array=pdos_temp
       end if
    end if


    return


  contains

    !===============================================================================
    subroutine pdos_analyse_substring(ctemp,species_proj)
      !===============================================================================
      ! This is a mindbendingly horrific exercise in book-keeping
      !===============================================================================
      use od_cell, only : num_species,atoms_species_num
      use od_io, only : maxlen, io_error
      implicit none


      character(len=maxlen), intent(inout) :: ctemp
      integer, optional, intent(out) :: species_proj

      integer, save :: offset=0

      character(len=maxlen) :: ctemp2,c_am,m_string

      integer :: pos_l,pos_r,ia,iz,idiff,ic1,ic2,species,num_sites,num_am
      character(len=3)  :: c_symbol='   '
      logical :: am_sum,site_sum

      integer   :: num1,num2,i_punc,pos3,loop_l,loop_a,loop_p,loop_j
      integer   :: counter,loop_r,range_size,ierr
      character(len=maxlen) :: dummy
      character(len=10), parameter :: c_digit="0123456789"
      character(len=1) , parameter :: c_range="-"
      character(len=1) , parameter :: c_sep=","
      character(len=4) , parameter :: c_punc=" ,-:"
      character(len=5)  :: c_num1,c_num2
      integer, allocatable :: pdos_atoms(:),pdos_ang(:)
      logical :: lcount


      allocate(pdos_atoms(maxval(atoms_species_num)),stat=ierr)
      if(ierr/=0) call io_error('Error: pdos_analyse_substring - allocation of pdos_atoms failed')
      allocate(pdos_ang(max_am),stat=ierr)
      if(ierr/=0) call io_error('Error: pdos_analyse_substring - allocation of pdos_ang failed')

      lcount=.false.
      if(present(species_proj)) lcount=.true.

      pdos_atoms=0;pdos_ang=0
      site_sum=.false.
      am_sum=.false.
      ! look for Ang Mtm string eg (s,p)
      c_am=''
      pos_l=index(ctemp,'(')
      if(pos_l>0) then
         pos_r=index(ctemp,')')
         if(pos_r==0) call io_error ('pdos_analyse_substring: found ( but no )')
         if(pos_r<=pos_l) call io_error ('pdos_analyse_substring: found ) before (')
         c_am=ctemp(pos_l+1:pos_r-1)
         ctemp=ctemp(:pos_l-1)
      else ! implicit sum over AM
         am_sum=.true.
      end if

      ia = ichar('a')
      iz = ichar('z')
      idiff = ichar('Z')-ichar('z')

      ic1=ichar(ctemp(1:1))
      if(ic1<ia .or. ic1>iz) call io_error ('pdos_analyse_substring: problem reading atomic symbol in pdos string')
      ic2=ichar(ctemp(2:2))
      if(ic2>=ia .and. ic1<=iz) then
         c_symbol(1:1)=char(ic1+idiff)
         c_symbol(2:2)=ctemp(2:2)
         ctemp=ctemp(3:)
      else
         c_symbol(1:1)=char(ic1+idiff)
         ctemp=ctemp(2:)
      end if
      species=0
      do loop_j=1,num_species
         if(adjustl(c_symbol)==adjustl(pdos_symbol(loop_j))) then
            species=loop_j
         end if
      end do
      if(species==0) call io_error('pdos_analyse_substring: Failed to match atomic symbol in pdos string')


      !Count atoms numbers
      counter=0
      dummy=adjustl(ctemp)
      if (len_trim(dummy)>0) then
         dummy=adjustl(dummy)
         do 
            i_punc=scan(dummy,c_punc)
            if(i_punc==0) call io_error('pdos_analyse_substring: error looking for atom numbers') 
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
                  if(min(num1,num2)+loop_r-1>pdos_sites(species)) &
                       call io_error&
          ('pdos_analyse_substring: Atom number given in pdos string is greater than number of atoms for given species')
                  pdos_atoms(min(num1,num2)+loop_r-1)=1
               end do
            else
               counter=counter+1
               if(num1>pdos_sites(species)) &
                    call io_error&
          ('pdos_analyse_substring: Atom number given in pdos string is greater than number of atoms for given species')
               pdos_atoms(num1)=1
            end if

            if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
            if(scan(dummy,c_range)==1) call io_error('pdos_analyse_substring: Error parsing atoms numbers - incorrect range') 
            if(index(dummy,' ')==1) exit
         end do
      else
         site_sum=.true.
      end if

      ! count am
      counter=0
      dummy=adjustl(c_am)
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
               call io_error('pdos_analyse_substring: Problem reading l state ')
            end select
            if (pos3==0) exit
            dummy=dummy(pos3+1:)
         enddo
      else
         am_sum=.true.
      endif

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
      if(lcount) species_proj=num_am*num_sites

      if(.not. lcount) then
         loop_p=1+offset
         if(site_sum.and.am_sum) then
            pdos_projection_array(species,:,:,loop_p)=1
         elseif(site_sum.and..not.am_sum) then
            do loop_l=1,max_am
               if(pdos_ang(loop_l)==0) cycle 
               pdos_projection_array(species,:,loop_l,loop_p)=1
               loop_p=loop_p+1
            end do
         elseif(.not.site_sum.and.am_sum) then
            do loop_a=1,pdos_sites(species)
               if(pdos_atoms(loop_a)==0) cycle 
               pdos_projection_array(species,loop_a,:,loop_p)=1
               loop_p=loop_p+1
            end do
         else
            do loop_l=1,max_am
               if(pdos_ang(loop_l)==0) cycle 
               do loop_a=1,pdos_sites(species)
                  if(pdos_atoms(loop_a)==0) cycle 
                  pdos_projection_array(species,loop_a,loop_l,loop_p)=1
                  loop_p=loop_p+1
               end do
            end do
         end if
         offset=loop_p-1
      end if


      return

101   call io_error('pdos_analyse_substring Error parsing keyword ')
106   call io_error('pdos_analyse_substring: Problem reading l state into string ')


    end subroutine pdos_analyse_substring


  end subroutine pdos_get_string


  subroutine pdos_analyse_orbitals
    use od_electronic, only :  pdos_orbital,pdos_mwab
    use od_cell, only : atoms_symbol,num_species
    use od_constants, only : periodic_table_name
    use od_io, only : io_error
    implicit none

    integer :: loop,loop2,counter,ierr

    if(maxval(pdos_orbital%species_no(:))>num_species) &
         call io_error('Error: pdos_analyse_substring - more species in pdos file than in cell file')

    allocate(pdos_sites(maxval(pdos_orbital%species_no(:))),stat=ierr)
    if(ierr/=0) call io_error('Error: pdos_analyse_substring - allocation of pdos_sites failed')
    allocate(pdos_am(maxval(pdos_orbital%species_no(:)),max_am),stat=ierr)
    if(ierr/=0) call io_error('Error: pdos_analyse_substring - allocation of pdos_am failed')
    allocate(pdos_symbol(maxval(pdos_orbital%species_no(:))),stat=ierr)
    if(ierr/=0) call io_error('Error: pdos_analyse_substring - allocation of pdos_symbol failed')
    pdos_sites=0;pdos_am=0
    do loop=1,pdos_mwab%norbitals
       if(pdos_orbital%rank_in_species(loop)>pdos_sites(pdos_orbital%species_no(loop))) &
            pdos_sites(pdos_orbital%species_no(loop))=pdos_orbital%rank_in_species(loop)
       if(pdos_orbital%rank_in_species(loop)==1) &
            pdos_am(pdos_orbital%species_no(loop),pdos_orbital%am_channel(loop)+1)=1
    end do

    !Now need to figure out symbols for each species

    counter=1
    do loop2=1,109
       do loop=1,num_species
          if(atoms_symbol(loop)==periodic_table_name(loop2)) then
             pdos_symbol(counter)=periodic_table_name(loop2)
             counter=counter+1
             !check atom count here
          end if
       end do
    end do


  end subroutine pdos_analyse_orbitals




  !===============================================================================
  subroutine pdos_write
    !===============================================================================
    ! Write out the pdos that was requested. Write them all to the same file, unless
    ! we don't have a short cut. In this case, write 10 projectors per file.
    !===============================================================================
    use od_io, only : seedname
    implicit none
   
    character(len=20) :: start_iproj_name, end_iproj_name 
    integer            ::  ifile, nfile, start_iproj, end_iproj 
    character(len=30)  :: name
    

   ! write(*,*) "======================================================================================"
   ! write(*,*) "ispecies,ispecies_num,iam,iproj,pdos_projection_array(ispecies,ispecies_num,iam,iproj)"
   ! do iproj=1,num_proj
   !    do iam=1,max_am
   !       do ispecies_num=1,maxval(atoms_species_num)
   !          do  ispecies=1,num_species   
   !             write(*,*) ispecies,ispecies_num,iam,iproj,pdos_projection_array(ispecies,ispecies_num,iam,iproj)
   !          enddo
   !       enddo
   !    enddo
   ! enddo
   ! write(*,*) "======================================================================================"

    if(shortcut) then
       ! write everything to one file
       name=trim(seedname)//'.pdos.dat'
       call write_proj_to_file(1,num_proj,name)
    else ! not shortcut
       nfile=int(num_proj/10)+1 ! Number of output files
       do ifile=1,nfile
          start_iproj=(ifile-1)*10+1 ! First projector in nfile
          if(ifile==nfile) then ! We're doing the last file
             end_iproj=num_proj
          else
             end_iproj=ifile*10
          endif
          
          write(start_iproj_name,'(I20.4)') start_iproj
          write(end_iproj_name,'(I20.4)') end_iproj
          name=trim(seedname)//'.pdos.proj-'//trim(adjustl(start_iproj_name))//'-'//trim(adjustl(end_iproj_name))//'.dat'
          call write_proj_to_file(start_iproj,end_iproj,name) 
        enddo
     endif
   end subroutine pdos_write
   
   
   

  subroutine write_proj_to_file(start_proj,stop_proj,name)
    !===============================================================================
    ! Write out projectors, start_proj, stop_proj, to file name
    !===============================================================================
    use od_dos_utils,       only : E
    use od_parameters,only : dos_nbins, iprint
    use od_algorithms, only : channel_to_am
    use od_electronic, only         : pdos_mwab
    use od_cell, only : atoms_species_num, num_species 
    use od_io, only : io_file_unit, io_error, io_date, stdout
    
    
    implicit none
    integer, intent(in) :: start_proj, stop_proj
    character(len=30), intent(in) :: name
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: string 
    integer :: iproj, iam, ispecies_num, ispecies
    integer :: idos, i, pdos_file,ierr 
   
    write(string,'(I4,"(x,es14.7)")') (stop_proj-start_proj)+1

    pdos_file=io_file_unit()
    open(unit=pdos_file,file=trim(name),iostat=ierr)
    if(iprint>2) write(stdout,'(1x,30a,30a)') " Writing PDOS projectors to: ", trim(name)
    if(ierr.ne.0) call io_error(" ERROR: Cannot open output file in pdos: pdos_write")
    
    write(pdos_file, *) "##############################################################################"
    write(pdos_file,*) "#"
    write(pdos_file, *) "#                  O p t a D O S   o u t p u t   f i l e "  
    write(pdos_file, '(1x,a1)') "#"
    call io_date(cdate,ctime)
    write(pdos_file,*)  '#  Generated on ',cdate,' at ',ctime
    write(pdos_file, '(1x,a78)') "##############################################################################"
    write(pdos_file,'(1a,a)') '#','*----------------------------------------------------------------------------*'   
    write(pdos_file,'(1a,a)') '#','|                    Partial Density of States -- Projectors                 |'
    write(pdos_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'
    
    
    if(pdos_mwab%nspins>1) then
       do iproj=start_proj,stop_proj
          write(pdos_file,'(1a,a1,a12,i4,a10,50x,a1)') '#','|', ' Column: ',iproj, ' contains:', '|'
          write(pdos_file,'(1a,a1,a16,10x,a14,5x,a15,16x,a1)') '#','|', ' Atom ', ' AngM Channel ', ' Spin Channel ', '|'
             do  ispecies=1,num_species   
                do ispecies_num=1,atoms_species_num(ispecies)
                   do iam=1,max_am
                   if(pdos_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                      write(pdos_file,'(1a,a1,a13,i3,a18,16x,a2,24x,1a)') "#","|", pdos_symbol(ispecies), &
                           &ispecies_num, channel_to_am(iam),'Up','|'
                   endif
                enddo
             enddo
          enddo
          write(pdos_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'
       enddo
       do iproj=start_proj,stop_proj
          write(pdos_file,'(1a,a1,a12,i4,a10,50x,a1)') '#','|', ' Column: ',iproj+num_proj, ' contains:', '|'
          write(pdos_file,'(1a,a1,a16,10x,a14,5x,a15,16x,a1)') '#','|', ' Atom ', ' AngM Channel ', ' Spin Channel ', '|'
          do  ispecies=1,num_species   
             do ispecies_num=1,atoms_species_num(ispecies)
                do iam=1,max_am
                   if(pdos_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                      write(pdos_file,'(1a,a1,a13,i3,a18,15x,a4,23x,1a)') "#","|", pdos_symbol(ispecies), &
                           &ispecies_num, channel_to_am(iam),'Down','|'
                   endif
                enddo
             enddo
          enddo
          write(pdos_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'
       enddo
       
       dos_partial(:,2,:)=-dos_partial(:,2,:)
       do idos = 1,dos_nbins
          write(pdos_file,'(es14.7,'//trim(string)//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=start_proj,stop_proj) &
               & ,(dos_partial(idos,2,i),i=start_proj,stop_proj)
       end do
    else
       do iproj=start_proj,stop_proj
          write(pdos_file,'(1a,a1,a12,i4,a10,50x,a1)') '#','|', ' Projector: ',iproj, ' contains:', '|'
          write(pdos_file,'(1a,a1,a16,10x,a14,36x,a1)') '#','|', ' Atom ', ' AngM Channel ', '|'
          do  ispecies=1,num_species   
             do ispecies_num=1,atoms_species_num(ispecies)
                do iam=1,max_am
                   if(pdos_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                      write(pdos_file,'(1a,a1,a13,i3,a18,42x,a1)') "#","|", pdos_symbol(ispecies), &
                           &ispecies_num, channel_to_am(iam),'|' 
                   endif
                enddo
             enddo
          enddo
          write(pdos_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'
       enddo
       
       do idos = 1,dos_nbins
          write(pdos_file,'(es14.7,'//trim(string)//')') E(idos),(dos_partial(idos,1,i),i=start_proj,stop_proj)
       end do
    endif
    
    close(pdos_file)
  end subroutine write_proj_to_file
  
  subroutine pdos_report_projectors
    use od_algorithms, only : channel_to_am
    use od_cell, only : atoms_species_num, num_species 
    use od_io, only : stdout
    implicit none
    
    integer :: iproj, iam, ispecies_num, ispecies

    write(stdout,*)
    write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'   
    write(stdout,'(1x,a)') '|                    Partial Density of States -- Projectors                 |'
    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
    do iproj=1,num_proj
       write(stdout,'(1x,a1,a12,i4,a10,50x,a1)') '|', ' Projector: ',iproj, ' contains:', '|'
       write(stdout,'(1x,a1,a16,10x,a14,36x,a1)') '|', ' Atom ', ' AngM Channel ', '|'
     
       do  ispecies=1,num_species 
          do ispecies_num=1,atoms_species_num(ispecies)
             do iam=1,max_am
                if(pdos_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                   write(stdout,'(1x,a1,a13,i3,a18,42x,a1)') "|", pdos_symbol(ispecies), &
                        ispecies_num, channel_to_am(iam),'|' !, " |  DEBUG :",  ispecies ,iam
                endif
             enddo
          enddo
       enddo
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
    enddo
  end subroutine pdos_report_projectors
  
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


!!$  subroutine general_write_pdos
!!$    !===============================================================================
!!$    ! Write out the pdos that was requested. Make a pretty header so that the user
!!$    ! knows what each column means
!!$    !===============================================================================
!!$    use od_dos_utils,       only : E
!!$    use od_parameters,only : dos_nbins
!!$    use od_algorithms, only : channel_to_am
!!$    use od_electronic, only         : pdos_mwab
!!$    use od_cell, only : atoms_species_num, num_species 
!!$    use od_io, only : io_file_unit, seedname, io_error, io_date
!!$
!!$    implicit none
!!$
!!$   character(len=11) :: cdate
!!$   character(len=9) :: ctime
!!$    character(len=20) :: string, filename 
!!$    integer :: iproj, iam, ispecies_num, ispecies, species, species_num
!!$    integer :: last_species, last_species_num
!!$    integer :: idos, i, pdos_file,ierr, start_proj 
!!$
!!$    logical :: projector_to_file
!!$
!!$
!!$    write(string,'(I4,"(x,es14.7)")') pdos_mwab%norbitals
!!$
!!$    start_proj=1
!!$    projectors: do iproj=1,num_proj
!!$       projector_to_file=.false.
!!$     
!!$       ! Are we writing .pdos.projX.dat or .pdos.AtomAtomNo.dat?
!!$       ! does this projector contain more than one atom?
!!$       do iam=1,max_am       
!!$          if(sum(pdos_projection_array(:,:,iam,iproj))>1) then
!!$             ! Yes it does contain more than one atom
!!$             projector_to_file=.true.
!!$          endif
!!$       enddo
!!$
!!$       if(projector_to_file) then
!!$          ! Then let's write out this projector and move on to the next one
!!$
!!$          ! Must first check whether this isn't the last one in a previous projector group
!!$          if(start_proj.ne.iproj) then ! Yes it is.
!!$             write(string,'(I20)') last_species_num
!!$             filename=pdos_symbol(last_species)//adjustl(string)
!!$             write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ", &
!!$                  & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$             call write_proj_to_file(start_proj, iproj, filename)
!!$             start_proj=iproj+1 ! Reset start counter, and we've written the current one.
!!$             cycle projectors
!!$          endif
!!$
!!$          write(*,*) "For proj:", iproj, "there is more than one atom"
!!$          write(*,*) "Hence we're writing projectors to files"
!!$          write(string,'(I20)') iproj
!!$          filename="proj"//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ", &
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(iproj,iproj,filename) 
!!$          start_proj=iproj+1 ! Reset start counter, and we've written the current one.
!!$          cycle projectors
!!$       endif
!!$
!!$       ! Since we're not writing just one projector to the file, we're going to have to work out
!!$       ! how many projectors there are going to be in this file.
!!$       ! Work out the species and species_rank of this projector
!!$       scan: do iam=1,max_am
!!$          do ispecies_num=1,maxval(atoms_species_num)
!!$             do  ispecies=1,num_species   
!!$                if(pdos_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
!!$                   species=ispecies
!!$                   species_num=ispecies_num
!!$                   write(*,*) "Projector ",iproj," is Species ", ispecies, " Rank ", ispecies_num 
!!$                   exit scan
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo scan
!!$    
!!$       ! First time through we just put the info about this projector into the registry
!!$       if(iproj==start_proj) then
!!$          write(*,*) "Skipping over projector:", iproj," as we've nothing to compare it against yet"
!!$          last_species=species
!!$          last_species_num=species_num
!!$          start_proj=iproj
!!$       ! If this is the same species as the last one. We go around again.
!!$       elseif((species==last_species).and.(species_num==last_species_num)) then
!!$          write(*,*) "Projector ", iproj, " has the same Species and Rank as ", start_proj
!!$          last_species=species
!!$          last_species_num=species_num
!!$       else ! We've come to the end of the projector group, so need to write out all the old ones.
!!$          write(*,*) "Projector ", iproj, " has different Species and Rank to ", start_proj
!!$          write(string,'(I20)') iproj
!!$          filename=trim(pdos_symbol(species))//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj-1, " to file ", &
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(start_proj, iproj-1,filename)
!!$          start_proj=iproj ! Since we haven't written the current projector yet
!!$          last_species=species
!!$          last_species_num=species_num
!!$       endif
!!$
!!$       ! If this is our last loop, then we'd better write the last on out too.
!!$       if(iproj==num_proj) then
!!$          write(*,*) "Last Projector group ", start_proj, " to ", iproj
!!$          write(string,'(I20)') iproj
!!$          filename=trim(pdos_symbol(species))//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ",&
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(start_proj, iproj, filename)
!!$       endif
!!$
!!$    enddo projectors
!!$    return
!!$  end subroutine general_write_pdos


end module od_pdos
