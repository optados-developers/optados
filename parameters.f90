!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! This module contains GPL routines from Wannier90           !
! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          !
!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This version (c) Jonathan Yates 2010                       !
!                                                            !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module od_parameters

  use constants, only : dp
  use od_io,        only : stdout,maxlen

  implicit none

  private

  !Input
  integer,           public, save :: iprint
  character(len=20), public, save :: energy_unit
  character(len=20), public, save :: length_unit
  integer,           public, save :: mp_grid(3)
  real(kind=dp),     public, save :: dis_win_min
  real(kind=dp),     public, save :: dis_win_max
  real(kind=dp),     public, save :: dis_froz_min

  integer,           public, save :: dos_num_points
  real(kind=dp),     public, save :: dos_energy_step
  real(kind=dp),     public, save :: dos_gaussian_width
  character(len=20), public, save :: dos_plot_format



  real(kind=dp), allocatable,    public, save :: kpt_latt(:,:) !kpoints in lattice vecs
  real(kind=dp),     public, save :: real_lattice(3,3)

  ! Atom sites
  real(kind=dp), allocatable,     public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp), allocatable,     public, save :: atoms_pos_cart(:,:,:)
  integer, allocatable,           public, save :: atoms_species_num(:)  
  character(len=maxlen), allocatable,  public, save :: atoms_label(:)
  character(len=2), allocatable,  public, save :: atoms_symbol(:)
  integer,                        public, save :: num_atoms
  integer,                        public, save :: num_species

  !parameters dervied from input
  integer,           public, save :: num_kpts
  real(kind=dp),     public, save :: recip_lattice(3,3)
  real(kind=dp),     public, save :: cell_volume
  real(kind=dp),     public, save :: real_metric(3,3)
  real(kind=dp),     public, save :: recip_metric(3,3)
  integer,           public, save :: bands_num_spec_points  
  character(len=1), allocatable,    public, save ::bands_label(:)
  real(kind=dp), allocatable,    public, save ::bands_spec_points(:,:)
  real(kind=dp), allocatable,    public, save ::kpt_cart(:,:) !kpoints in cartesians
  real(kind=dp),     public, save :: lenconfac

  integer :: num_lines
  character(len=maxlen), allocatable :: in_data(:)


  public :: param_read
  public :: param_write
  public :: param_dealloc


contains



  !==================================================================!
  subroutine param_read ( )
  !==================================================================!
  !                                                                  !
  ! Read parameters and calculate derived values                     !
  !                                                                  !
  !===================================================================  
    use constants, only : bohr
!    use w90_utility,   only : utility_recip_lattice,utility_compute_metric
    use od_io,        only : io_error,io_file_unit,seedname
    implicit none

    !local variables
    real(kind=dp)  :: real_lattice_tmp(3,3)
    integer :: nkp,i,j,n,k,itmp,i_temp,i_temp2,eig_unit,loop,ierr,iv_temp(3)
    logical :: found,found2,eig_found,lunits,chk_found
    character(len=6) :: spin_str
    real(kind=dp) :: cosa(3),rv_temp(3)

    call param_in_file

    iprint          =  1             ! Verbosity
    call param_get_keyword('iprint',found,i_value=iprint)

    energy_unit     =  'ev'          !
    call param_get_keyword('energy_unit',found,c_value=energy_unit)

    length_unit     =  'ang'         !
    lenconfac=1.0_dp
    call param_get_keyword('length_unit',found,c_value=length_unit)
    if (length_unit.ne.'ang' .and. length_unit.ne.'bohr') &
         call io_error('Error: value of length_unit not recognised in param_read')
    if (length_unit.eq.'bohr') lenconfac=1.0_dp/bohr

    dos_num_points            = 50
    call param_get_keyword('dos_num_points',found,i_value=dos_num_points)
    if (dos_num_points<0) call io_error('Error: dos_num_points must be positive')       

    dos_energy_step           = 0.01_dp
    call param_get_keyword('dos_energy_step',found,r_value=dos_energy_step)

    dos_gaussian_width        = 0.1_dp
    call param_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

    dos_plot_format           = 'gnuplot'
    call param_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    call param_uppercase()

303  continue

    deallocate(in_data,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating in_data in param_read')

    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !


    return

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')
106 call io_error('Error: Problem reading eigenvalue file '//trim(seedname)//'.eig')

  end subroutine param_read


  !===================================================================
  subroutine param_uppercase
    !===================================================================
    !                                                                  !
    ! Convert a few things to uppercase to look nice in the output     !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: nsp,ic,loop

    ! Atom labels (eg, si --> Si)
    do nsp=1,num_species
       ic=ichar(atoms_label(nsp)(1:1))                           
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_label(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
    enddo

    do nsp=1,num_species
       ic=ichar(atoms_symbol(nsp)(1:1))
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_symbol(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
    enddo


    ! Bands labels (eg, x --> X)
    do loop=1,bands_num_spec_points
       ic=ichar(bands_label(loop))                           
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            bands_label(loop) = char(ic+ichar('Z')-ichar('z'))
    enddo

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic=ichar(length_unit(1:1))
    if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
         length_unit(1:1) = char(ic+ichar('Z')-ichar('z'))

    return

  end subroutine param_uppercase


  !===================================================================
  subroutine param_write
    !==================================================================!
    !                                                                  !
    ! write parameters to stdout                                       !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: i,nkp,loop,nat,nsp

    ! System
    write(stdout,*)
    write(stdout,'(36x,a6)') '------' 
    write(stdout,'(36x,a6)') 'SYSTEM' 
    write(stdout,'(36x,a6)') '------' 
    write(stdout,*)
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(30x,a21)') 'Lattice Vectors (Ang)' 
    else
       write(stdout,'(28x,a22)') 'Lattice Vectors (Bohr)' 
    endif
    write(stdout,101) 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    write(stdout,*)   
    write(stdout,'(19x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(2x,a7)') '(Ang^3)'
    else
       write(stdout,'(2x,a8)') '(Bohr^3)'
    endif
    write(stdout,*)   
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
       write(stdout,'(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write(stdout,101) 'b_1',(recip_lattice(1,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_2',(recip_lattice(2,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_3',(recip_lattice(3,I)/lenconfac, i=1,3)
    write(stdout,*)   ' '
    ! Atoms
    if(num_atoms>0) then
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
       else
          write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
       endif
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(stdout,'(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp),nat,atoms_pos_frac(:,nat,nsp),&
                  '|',atoms_pos_cart(:,nat,nsp)*lenconfac,'|'
          end do
       end do
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
    else
       write(stdout,'(25x,a)') 'No atom positions specified'
    end if
    write(stdout,*) ' '

    ! K-points
    write(stdout,'(32x,a)') '------------'
    write(stdout,'(32x,a)') 'K-POINT GRID'
    write(stdout,'(32x,a)') '------------'
    write(stdout,*) ' '
    write(stdout,'(13x,a,i3,1x,a1,i3,1x,a1,i3,6x,a,i5)') 'Grid size =',mp_grid(1),'x',mp_grid(2),'x',mp_grid(3),&
         'Total points =',num_kpts
    write(stdout,*) ' '
    if(iprint>1) then
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
       else
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
       endif
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nkp=1,num_kpts
          write(stdout,'(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',nkp,kpt_latt(:,nkp),'|',kpt_cart(:,nkp)/lenconfac,'|'
       end do
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       write(stdout,*) ' '
    end if
    ! Main
    write(stdout,*) ' '
    write(stdout,'(1x,a78)') '*---------------------------------- MAIN ------------------------------------*'
    write(stdout,'(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :',trim(length_unit),'|'  
    write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'

    !
    write(stdout,'(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
    !


101 format(20x,a3,2x,3F11.6)

  end subroutine param_write


  !==================================================================!
  subroutine param_dealloc
    !==================================================================!
    !                                                                  !
    ! release memory from allocated parameters                         !
    !                                                                  !
    !===================================================================  
    use od_io, only : io_error

    implicit none
    integer :: ierr

!    if ( allocated ( ndimwin ) ) then
!       deallocate (  ndimwin, stat=ierr  )
!       if (ierr/=0) call io_error('Error in deallocating ndimwin in param_dealloc')
!    end if
    return

  end subroutine param_dealloc



  !=======================================!
  subroutine param_in_file
    !=======================================!
    ! Load the *.win file into a character  !
    ! array in_file, ignoring comments and  !
    ! blank lines and converting everything !
    ! to lowercase characters               !
    !=======================================!

    use od_io,        only : io_file_unit,io_error,seedname
    use algorithms,   only : utility_lowercase

    implicit none

    integer           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy

    in_unit=io_file_unit( )
    open (in_unit, file=trim(seedname)//'.win',form='formatted',status='old',err=101)

    num_lines=0;tot_num_lines=0
    do
       read(in_unit, '(a)', iostat = ierr, err= 200, end =210 ) dummy
       dummy=adjustl(dummy)
       tot_num_lines=tot_num_lines+1
       if( .not.dummy(1:1)=='!'  .and. .not. dummy(1:1)=='#' ) then
          if(len(trim(dummy)) > 0 ) num_lines=num_lines+1
       endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.win')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.win')
210 continue
    rewind(in_unit)

    allocate(in_data(num_lines),stat=ierr)
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

    line_counter=0
    do loop=1,tot_num_lines
       read(in_unit, '(a)', iostat = ierr, err= 200 ) dummy
       dummy=utility_lowercase(dummy)
       dummy=adjustl(dummy)
       if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
       if(len(trim(dummy)) == 0 ) cycle
       line_counter=line_counter+1
       in1=index(dummy,'!')
       in2=index(dummy,'#')
       if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
       if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
       if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
       if(in2>0 .and. in1>0 )   in_data(line_counter)=dummy(:min(in1,in2)-1)
    end do

    close(in_unit)

  end subroutine param_in_file


  !===========================================================================!
  subroutine param_get_keyword(keyword,found,c_value,l_value,i_value,r_value)
    !===========================================================================!
    !                                                                           !
    !             Finds the value of the required keyword.                      !
    !                                                                           !
    !===========================================================================!

    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    character(*)     ,optional, intent(inout) :: c_value
    logical          ,optional, intent(inout) :: l_value
    integer          ,optional, intent(inout) :: i_value
    real(kind=dp)    ,optional, intent(inout) :: r_value

    integer           :: kl, in,loop,itmp
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.

    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       itmp=in+len(trim(keyword))
       if (in_data(loop)(itmp:itmp)/='=' &
            .and. in_data(loop)(itmp:itmp)/=':' &
            .and. in_data(loop)(itmp:itmp)/=' ') cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(found) then
       if( present(c_value) ) c_value=dummy
       if( present(l_value) ) then
          if (index(dummy,'t') > 0) then
             l_value=.true.
          elseif (index(dummy,'f') > 0) then
             l_value=.false.
          else
             call io_error('Error: Problem reading logical keyword '//trim(keyword))
          endif
       endif
       if( present(i_value) ) read(dummy,*,err=220,end=220) i_value
       if( present(r_value) ) read(dummy,*,err=220,end=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))


  end subroutine param_get_keyword


  !=========================================================================================!
  subroutine param_get_keyword_vector(keyword,found,length,c_value,l_value,i_value,r_value)
    !=========================================================================================!
    !                                                                                         !
    !                  Finds the values of the required keyword vector                        !
    !                                                                                         !
    !=========================================================================================!

    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(in)  :: length
    character(*)     ,optional, intent(inout) :: c_value(length)
    logical          ,optional, intent(inout) :: l_value(length)
    integer          ,optional, intent(inout) :: i_value(length)
    real(kind=dp)    ,optional, intent(inout) :: r_value(length)

    integer           :: kl, in,loop,i
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.



    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(found) then
       if( present(c_value) ) read(dummy,*,err=230,end=230) (c_value(i),i=1,length)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=230,end=230) (i_value(i),i=1,length)
       if( present(r_value) ) read(dummy,*,err=230,end=230) (r_value(i),i=1,length)
    end if



    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector')


  end subroutine param_get_keyword_vector



  !========================================================!
  subroutine param_get_vector_length(keyword,found,length)
    !======================================================!
    !                                                      !
    !        Returns the length of a keyword vector        !
    !                                                      !
    !======================================================!

    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(out)  :: length

    integer           :: kl, in,loop,pos
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.



    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    length=0
    if(found) then
       if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
       length=1
       dummy=adjustl(dummy)
       do 
          pos=index(dummy,' ')
          dummy=dummy(pos+1:)
          dummy=adjustl(dummy)
          if(len_trim(dummy)>0) then
             length=length+1
          else
             exit
          endif

       end do

    end if



    return


  end subroutine param_get_vector_length


  !==============================================================================================!
  subroutine param_get_keyword_block(keyword,found,rows,columns,c_value,l_value,i_value,r_value)
    !==============================================================================================!
    !                                                                                              !
    !                           Finds the values of the required data block                        !
    !                                                                                              !
    !==============================================================================================!

    use constants, only : bohr
    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(in)  :: rows
    integer,           intent(in)  :: columns
    character(*)     ,optional, intent(inout) :: c_value(columns,rows)
    logical          ,optional, intent(inout) :: l_value(columns,rows)
    integer          ,optional, intent(inout) :: i_value(columns,rows)
    real(kind=dp)    ,optional, intent(inout) :: r_value(columns,rows)

    integer           :: in,ins,ine,loop,i,line_e,line_s,counter,blen
    logical           :: found_e,found_s,lconvert
    character(len=maxlen) :: dummy,end_st,start_st

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do

    if(.not. found_s) then
       found=.false.
       return
    end if


    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    ! number of lines of data in block
    blen = line_e-line_s-1

    !    if( blen /= rows) then
    !       if ( index(trim(keyword),'unit_cell_cart').ne.0 ) then
    !          if ( blen /= rows+1 ) call io_error('Error: Wrong number of lines in block '//trim(keyword))
    !       else
    !          call io_error('Error: Wrong number of lines in block '//trim(keyword))          
    !       endif
    !    endif

    if ( (blen.ne.rows) .and. (blen.ne.rows+1) ) &
         call io_error('Error: Wrong number of lines in block '//trim(keyword))          

    if ( (blen.eq.rows+1) .and. (index(trim(keyword),'unit_cell_cart').eq.0) ) &
         call io_error('Error: Wrong number of lines in block '//trim(keyword))          


    found=.true.

    lconvert=.false.
    if (blen==rows+1) then
       dummy=in_data(line_s+1)
       if ( index(dummy,'ang').ne.0 ) then
          lconvert=.false.
       elseif ( index(dummy,'bohr').ne.0 ) then
          lconvert=.true.
       else
          call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       endif
       in_data(line_s)(1:maxlen) = ' '
       line_s=line_s+1
    endif

!    r_value=1.0_dp
    counter=0
    do loop=line_s+1,line_e-1
       dummy=in_data(loop)
       counter=counter+1
       if( present(c_value) ) read(dummy,*,err=240,end=240) (c_value(i,counter),i=1,columns)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=240,end=240) (i_value(i,counter),i=1,columns)
       if( present(r_value) ) read(dummy,*,err=240,end=240) (r_value(i,counter),i=1,columns)
    end do

    if (lconvert) then
       if (present(r_value)) then
          r_value=r_value*bohr
       endif
    endif

    in_data(line_s:line_e)(1:maxlen) = ' '


    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))


  end subroutine param_get_keyword_block

  !=====================================================!
  subroutine param_get_block_length(keyword,found,rows,lunits)
    !=====================================================!
    !                                                     !
    !       Finds the length of the data block            !
    !                                                     !
    !=====================================================!

    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(out) :: rows
    logical, optional, intent(out) :: lunits

    integer           :: i,in,ins,ine,loop,line_e,line_s
    logical           :: found_e,found_s
    character(len=maxlen) :: end_st,start_st,dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)

    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do


    if(.not. found_s) then
       found=.false.
       return
    end if


    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do


    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    rows=line_e-line_s-1

    found=.true.

    if (present(lunits)) then
       dummy=in_data(line_s+1)
       !       write(stdout,*) dummy
       !       write(stdout,*) trim(dummy)
       read(dummy,*,end=555) atsym, (atpos(i),i=1,3)
       lunits=.false.
    endif

    if(rows<=0) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if


    return

555 lunits=.true.

    if(rows<=1) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if


    return

  end subroutine param_get_block_length




    !====================================================================!
    subroutine param_get_range_vector(keyword,found,length,lcount,i_value)
    !====================================================================!
    !   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100               !
    !   if(lcount) we return the number of states in length              !
    !====================================================================!
    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)    :: keyword
    logical          , intent(out)   :: found
    integer,           intent(inout) :: length
    logical,           intent(in)    :: lcount
    integer, optional, intent(out)   :: i_value(length)

    integer   :: kl, in,loop,num1,num2,i_punc
    integer   :: counter,i_digit,loop_r,range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit="0123456789"
    character(len=2) , parameter :: c_range="-:"
    character(len=3) , parameter :: c_sep=" ,;"
    character(len=5) , parameter :: c_punc=" ,;-:"
    character(len=5)  :: c_num1,c_num2

    
    if(lcount .and. present(i_value) ) call io_error('param_get_range_vector: incorrect call')

    kl=len_trim(keyword)

    found=.false.

    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       dummy=adjustl(dummy)
       if(.not. lcount) in_data(loop)(1:maxlen) = ' '
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(.not. found) return

    counter=0
    if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
    dummy=adjustl(dummy)
    do 
       i_punc=scan(dummy,c_punc)
       if(i_punc==0) call io_error('Error parsing keyword '//trim(keyword)) 
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
             if(.not. lcount) i_value(counter)=min(num1,num2)+loop_r-1
          end do
       else
          counter=counter+1 
          if(.not. lcount) i_value(counter)=num1
       end if

       if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
       if(scan(dummy,c_range)==1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range') 
       if(index(dummy,' ')==1) exit
    end do

    if(lcount) length=counter
    if(.not.lcount) then
       do loop=1,counter-1
          do loop_r=loop+1,counter 
             if(i_value(loop)==i_value(loop_r)) &
                call io_error('Error parsing keyword '//trim(keyword)//' duplicate values')
          end do
        end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword))


   end  subroutine param_get_range_vector


  !===================================!
  subroutine param_get_keyword_kpath
    !===================================!
    !                                   !
    !  Fills the kpath data block       !
    !                                   !
    !===================================!
    use od_io,        only : io_error

    implicit none

    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,i,line_e,line_s,counter
    logical           :: found_e,found_s
    character(len=maxlen) :: dummy,end_st,start_st

    keyword="kpoint_path"

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do



    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    counter=0
    do loop=line_s+1,line_e-1

       counter=counter+2
       dummy=in_data(loop)
       read(dummy,*,err=240,end=240) bands_label(counter-1),(bands_spec_points(i,counter-1),i=1,3)&
            ,bands_label(counter),(bands_spec_points(i,counter),i=1,3)
    end do


    in_data(line_s:line_e)(1:maxlen) = ' '

    return


240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end subroutine param_get_keyword_kpath



end module od_parameters
