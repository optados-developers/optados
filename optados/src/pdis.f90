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
module od_pdis
  !===============================================================================
  ! By exploiting the matrix weights read when performing an OptaDOS PDOS 
  ! calculation, this module implements the "PDIS" task, or projected dispersion/
  ! projected bandstructures. This is achieved in 3 steps:
  ! (1) Read the "pdos_weights/pdos_bin" file from a CASTEP spectral calculation
  !      on a kpoint path (as opposed to grid),
  ! (2) The exercise in book-keeping necessary to merge and discard the weights 
  !     such that the requested channels are given to the dos module,
  ! (3) The writing out of the pdis, one kpoint as a time, labelling the columns 
  !     for easy reading.
  !===============================================================================

  use od_constants, only : dp

  implicit none


  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  real(kind=dp), public, allocatable, save  :: dos_partial(:,:,:)
  real(kind=dp),allocatable, public, dimension(:,:,:,:) :: matrix_weights
  integer, allocatable :: pdis_projection_array(:,:,:,:)
  integer :: num_proj

  integer, parameter :: max_am=4                   ! s,p,d,f  hard coded!

  ! Data derived from the info in the pdos_weights file 
  character(len=3), allocatable :: pdos_symbol(:)   ! symbols
  integer, allocatable :: pdos_am(:,:)              ! angular mtm (num_species,max_am)
  integer, allocatable :: pdos_sites(:)             ! number of each species
  logical :: shortcut


  !-------------------------------------------------------------------------!

  private

  public :: pdis_calculate

  contains

  subroutine pdis_calculate
    use od_electronic, only : elec_pdos_read,efermi_set
    use od_dos_utils, only : dos_utils_calculate,dos_utils_set_efermi
    use od_comms, only : on_root
    use od_parameters, only : set_efermi_zero
    use od_io, only : stdout

    implicit none

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') '+                 Projected Dispersion Curve Calculation                     +'
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,'(1x,a78)') 
    endif

    ! read in the pdos weights
    call elec_pdos_read

    ! look at the orbitals and figure out which atoms / states we have
    call pdis_analyse_orbitals

    ! parse the pdis string to see what we want
    call pdis_get_string

    ! form the right matrix elements
    call pdis_merge

    ! and write everything out
    if(set_efermi_zero .and. .not.efermi_set) call dos_utils_set_efermi
    if (on_root) then
       call pdis_write
    endif


  end subroutine pdis_calculate

  !===============================================================================
  subroutine pdis_merge
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
    if(ierr/=0) call io_error('Error: pdis_merge - allocation of matrix_weights failed')
    matrix_weights=0.0_dp

    do N=1,num_kpoints_on_node(my_node_id) ! Loop over kpoints
       do N_spin=1,nspins                  ! Loop over spins
          do n_eigen=1,pdos_mwab%nbands    ! Loop over unoccupied states
             do nproj=1,num_proj
                do orb=1,pdos_mwab%norbitals
                   if(pdis_projection_array(pdos_orbital%species_no(orb),pdos_orbital%rank_in_species(orb) &
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

  end subroutine pdis_merge



  !===============================================================================
  subroutine pdis_get_string
    !===============================================================================
    ! This is a mindbendingly horrific exercise in book-keeping
    !===============================================================================
    use od_parameters, only : pdis_string
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

    ctemp=pdis_string
    if(index(ctemp,'species_ang')>0) then
       num_proj=0
       do loop=1,num_species
          num_proj=num_proj+count(pdos_am(loop,:)==1)
       end do
       allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')

       pdis_projection_array=0
       loop_p=1
       do loop=1,num_species
          do loop_l=1,max_am
             if(pdos_am(loop,loop_l)==0) cycle 
             pdis_projection_array(loop,:,loop_l,loop_p)=1
             loop_p=loop_p+1
          end do
       end do
       shortcut=.true.
    elseif(index(ctemp,'species')>0) then
       num_proj=num_species
       allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')

       pdis_projection_array=0
       do loop=1,num_species
          pdis_projection_array(loop,:,:,loop)=1
       end do
       shortcut=.true.
    elseif(index(ctemp,'sites')>0) then
       num_proj=0
       do loop=1,num_species
          num_proj=num_proj+pdos_sites(loop)
       end do
       allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')

       pdis_projection_array=0
       loop_p=1
       do loop=1,num_species
          do loop_a=1,pdos_sites(loop)
             pdis_projection_array(loop,loop_a,:,loop_p)=1
             loop_p=loop_p+1
          end do
       end do
       shortcut=.true.
    elseif(index(ctemp,'angular')>0) then
       num_proj=max_am
       allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')
       pdis_projection_array=0
       loop_p=0
       do loop=1,num_proj
          pdis_projection_array(:,:,loop,loop)=1
       end do
       shortcut=.true.
    endif


    if(.not. shortcut) then

       !look for sum
       ctemp=pdis_string
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
          call pdis_analyse_substring(ctemp3,species_proj)
          num_proj=num_proj+species_proj
          if (pos==0) exit
          species_count=species_count+1
          ctemp2=ctemp2(pos+1:)
       end do

       ! now allocate the correct sized array

       allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')
       pdis_projection_array=0

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
          call pdis_analyse_substring(ctemp3)
          ctemp2=ctemp2(pos+1:)
       end do

       if(pdos_sum) then
          allocate(pdos_temp(num_species,maxval(atoms_species_num),max_am,1),stat=ierr)
       if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdos_temp failed')

          pdos_temp=0
          do loop4=1,num_proj
             do loop3=1,max_am
                do loop2=1,maxval(atoms_species_num)
                   do loop=1,num_species
                      if(pdis_projection_array(loop,loop2,loop3,loop4)==1) then
                         pdos_temp(loop,loop2,loop3,1)=1
                      end if
                   end do
                end do
             end do
          end do
          deallocate(pdis_projection_array,stat=ierr)
          if(ierr/=0) call io_error('Error: pdis_get_string - deallocation of pdis_projection_array failed')
          num_proj=1
          allocate(pdis_projection_array(num_species,maxval(atoms_species_num),max_am,num_proj),stat=ierr)
          if(ierr/=0) call io_error('Error: pdis_get_string - allocation of pdis_projection_array failed')
          pdis_projection_array=0
          pdis_projection_array=pdos_temp
       end if
    end if


    return


  contains

    !===============================================================================
    subroutine pdis_analyse_substring(ctemp,species_proj)
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
      if(ierr/=0) call io_error('Error: pdis_analyse_substring - allocation of pdos_atoms failed')
      allocate(pdos_ang(max_am),stat=ierr)
      if(ierr/=0) call io_error('Error: pdis_analyse_substring - allocation of pdos_ang failed')

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
         if(pos_r==0) call io_error ('pdis_analyse_substring: found ( but no )')
         if(pos_r<=pos_l) call io_error ('pdis_analyse_substring: found ) before (')
         c_am=ctemp(pos_l+1:pos_r-1)
         ctemp=ctemp(:pos_l-1)
      else ! implicit sum over AM
         am_sum=.true.
      end if

      ia = ichar('a')
      iz = ichar('z')
      idiff = ichar('Z')-ichar('z')

      ic1=ichar(ctemp(1:1))
      if(ic1<ia .or. ic1>iz) call io_error ('pdis_analyse_substring: problem reading atomic symbol in pdos string')
      ic2=ichar(ctemp(2:2))
      if(ic2>=ia .and. ic1<=iz) then
         c_symbol(1:1)=char(ic1+idiff)
         c_symbol(2:2)=ctemp(2:2)
         ctemp=ctemp(3:)
      else
         c_symbol(1:1)=char(ic1+idiff)
         c_symbol(2:2)=''
         ctemp=ctemp(2:)
      end if
      species=0
      do loop_j=1,num_species
         if(adjustl(c_symbol)==adjustl(pdos_symbol(loop_j))) then
            species=loop_j
         end if
      end do
      if(species==0) call io_error('pdis_analyse_substring: Failed to match atomic symbol in pdos string')


      !Count atoms numbers
      counter=0
      dummy=adjustl(ctemp)
      if (len_trim(dummy)>0) then
         dummy=adjustl(dummy)
         do 
            i_punc=scan(dummy,c_punc)
            if(i_punc==0) call io_error('pdis_analyse_substring: error looking for atom numbers') 
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
          ('pdis_analyse_substring: Atom number given in pdos string is greater than number of atoms for given species')
                  pdos_atoms(min(num1,num2)+loop_r-1)=1
               end do
            else
               counter=counter+1
               if(num1>pdos_sites(species)) &
                    call io_error&
          ('pdis_analyse_substring: Atom number given in pdos string is greater than number of atoms for given species')
               pdos_atoms(num1)=1
            end if

            if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
            if(scan(dummy,c_range)==1) call io_error('pdis_analyse_substring: Error parsing atoms numbers - incorrect range') 
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
               call io_error('pdis_analyse_substring: Problem reading l state ')
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
            pdis_projection_array(species,:,:,loop_p)=1
         elseif(site_sum.and..not.am_sum) then
            do loop_l=1,max_am
               if(pdos_ang(loop_l)==0) cycle 
               pdis_projection_array(species,:,loop_l,loop_p)=1
               loop_p=loop_p+1
            end do
         elseif(.not.site_sum.and.am_sum) then
            do loop_a=1,pdos_sites(species)
               if(pdos_atoms(loop_a)==0) cycle 
               pdis_projection_array(species,loop_a,:,loop_p)=1
               loop_p=loop_p+1
            end do
         else
            do loop_l=1,max_am
               if(pdos_ang(loop_l)==0) cycle 
               do loop_a=1,pdos_sites(species)
                  if(pdos_atoms(loop_a)==0) cycle 
                  pdis_projection_array(species,loop_a,loop_l,loop_p)=1
                  loop_p=loop_p+1
               end do
            end do
         end if
         offset=loop_p-1
      end if


      return

101   call io_error('pdis_analyse_substring Error parsing keyword ')
106   call io_error('pdis_analyse_substring: Problem reading l state into string ')

    end subroutine pdis_analyse_substring

    

  end subroutine pdis_get_string


  subroutine pdis_analyse_orbitals
    use od_electronic, only :  pdos_orbital,pdos_mwab
    use od_cell, only : atoms_symbol,num_species
    use od_constants, only : periodic_table_name
    use od_io, only : io_error 
    implicit none

    integer :: loop,loop2,counter,ierr

    if(maxval(pdos_orbital%species_no(:))>num_species) then
         call io_error('Error: pdis_analyse_orbitals - more species in pdis file than in cell file')
    end if

    allocate(pdos_sites(maxval(pdos_orbital%species_no(:))),stat=ierr)
    if(ierr/=0) call io_error('Error: pdis_analyse_orbitals - allocation of pdis_sites failed')
    allocate(pdos_am(maxval(pdos_orbital%species_no(:)),max_am),stat=ierr)
    if(ierr/=0) call io_error('Error: pdis_analyse_orbitals - allocation of pdis_am failed')
    allocate(pdos_symbol(maxval(pdos_orbital%species_no(:))),stat=ierr)
    if(ierr/=0) call io_error('Error: pdis_analyse_orbitals - allocation of pdis_symbol failed')
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

  end subroutine pdis_analyse_orbitals
  
  !===============================================================================
  subroutine pdis_write
    !===============================================================================
    ! Write out the pdis that was requested. Write them all to the same file, kpoint
    ! by kpoint. 
    !===============================================================================
    use od_io, only : seedname
    implicit none
    character(len=512)  :: name
    
    ! write everything to one file
    name=trim(seedname)//'.pdis.dat'
    call write_pdis_to_file(1,num_proj,name)

   end subroutine pdis_write
   
  subroutine write_pdis_to_file(start_proj,stop_proj,name)
    !===============================================================================
    ! Write out projectors, start_proj, stop_proj, to file name, one kpoint at a time
    !===============================================================================
    use od_parameters, only : iprint
    use od_algorithms, only : channel_to_am
    use od_electronic, only : pdos_mwab, all_kpoints, band_energy
    use od_cell,       only : atoms_species_num, num_species, nkpoints
    use od_io,         only : io_file_unit, io_error, io_date, stdout
    
    implicit none
    integer, intent(in) :: start_proj, stop_proj
    character(len=512), intent(in) :: name
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: string 
    integer :: iproj, iam, ispecies_num, ispecies
    integer :: i, pdis_file,ierr 
    integer :: N, n_eigen

    write(string,'(I4,"(x,es14.7)")') (stop_proj-start_proj)+1

    pdis_file=io_file_unit()
    open(unit=pdis_file,file=trim(name),iostat=ierr)
    if(iprint>2) write(stdout,'(1x,a30,a30,17x,a1)') "| Writing PDIS projectors to: ", trim(name),"|"
    if(ierr.ne.0) call io_error(" ERROR: Cannot open output file in pdis: pdis_write")
    
    write(pdis_file, *) "##############################################################################"
    write(pdis_file,*) "#"
    write(pdis_file, *) "#                  O p t a D O S   o u t p u t   f i l e "  
    write(pdis_file, '(1x,a1)') "#"
    call io_date(cdate,ctime)
    write(pdis_file,*)  '#  Generated on ',cdate,' at ',ctime
    write(pdis_file, '(1x,a78)') "##############################################################################"
    write(pdis_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'   
    write(pdis_file,'(1a,a)') '#','|                    Projected Dispersion Curve -- Projectors                |'
    write(pdis_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'

    
    if(pdos_mwab%nspins>1) then
        call io_error('pDIS not implemented for multiple spin channels.')
    else
       do iproj=start_proj,stop_proj
          write(pdis_file,'(1a,a1,a12,i4,a10,50x,a1)') '#','|', ' Projector: ',iproj, ' contains:', '|'
          write(pdis_file,'(1a,a1,a16,10x,a14,36x,a1)') '#','|', ' Atom ', ' AngM Channel ', '|'
          do  ispecies=1,num_species   
             do ispecies_num=1,atoms_species_num(ispecies)
                do iam=1,max_am
                   if(pdis_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                      write(pdis_file,'(1a,a1,a13,i3,a18,42x,a1)') "#","|", pdos_symbol(ispecies), &
                           &ispecies_num, channel_to_am(iam),'|' 
                   endif
                enddo
             enddo
          enddo
          write(pdis_file,'(1a,a)') '#','+----------------------------------------------------------------------------+'
       enddo
       
       do N=1,nkpoints
          write(pdis_file, '(a10, i4, es14.7, es14.7, es14.7)') 'K-point   ', N, (all_kpoints(i, N), i=1,3)
          do n_eigen=1,pdos_mwab%nbands
             write(pdis_file, '(es20.7,'//trim(string)//')') band_energy(n_eigen, 1, N), &
                (matrix_weights(i, n_eigen, N, 1), i=start_proj,stop_proj)
          end do
        end do

    endif
    
    close(pdis_file)

  end subroutine write_pdis_to_file
  
  subroutine pdis_report_projectors
    use od_algorithms, only : channel_to_am
    use od_cell, only : atoms_species_num, num_species 
    use od_io, only : stdout
    implicit none
    
    integer :: iproj, iam, ispecies_num, ispecies

    write(stdout,*)
    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'   
    write(stdout,'(1x,a)') '|                    Projected Dispersion Curve -- Projectors                |'
    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
    do iproj=1,num_proj
       write(stdout,'(1x,a1,a12,i4,a10,50x,a1)') '|', ' Projector: ',iproj, ' contains:', '|'
       write(stdout,'(1x,a1,a16,10x,a14,36x,a1)') '|', ' Atom ', ' AngM Channel ', '|'
     
       do  ispecies=1,num_species 
          do ispecies_num=1,atoms_species_num(ispecies)
             do iam=1,max_am
                if(pdis_projection_array(ispecies,ispecies_num,iam,iproj)==1) then
                   write(stdout,'(1x,a1,a13,i3,a18,42x,a1)') "|", pdos_symbol(ispecies), &
                        ispecies_num, channel_to_am(iam),'|' !, " |  DEBUG :",  ispecies ,iam
                endif
             enddo
          enddo
       enddo
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
    enddo
  end subroutine pdis_report_projectors
  
end module od_pdis
