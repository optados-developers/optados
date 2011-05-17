!=========================================================================!
! Module: Cell                                                            !
! For dealing with the real and reciprocal space cell                     !
!-------------------------------------------------------------------------!
! Modules used:  constants                                                !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below                                                         !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris (So far)                            11/10/2010 !
!=========================================================================!
module od_cell
  use od_constants,only : dp
  use od_io,        only : maxlen 
  implicit none

  private
  !-------------------------------------------------------------------------!
  ! R E A L   S P A C E   V A R I A B L E S
  real(kind=dp), public, save :: real_lattice(1:3,1:3)
  real(kind=dp), public, save :: recip_lattice(1:3,1:3)
  real(kind=dp), public, save :: cell_volume
  !-------------------------------------------------------------------------!


  !-------------------------------------------------------------------------!
  ! R E C I P R O C A L   S P A C E   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: kpoint_r(:,:)
  real(kind=dp), allocatable, public, save :: kpoint_r_cart(:,:)
  real(kind=dp), allocatable, public, save :: kpoint_weight(:)
  integer, allocatable,public, save       :: num_kpoints_on_node(:)

  integer, public, save :: nkpoints 
  integer, public, save :: kpoint_grid_dim(3)
  !-------------------------------------------------------------------------!
  ! Symmetry Operations
  integer, public, save :: num_crystal_symmetry_operations
  real(kind=dp), allocatable, public, save :: crystal_symmetry_disps(:,:)
  real(kind=dp), allocatable, public, save :: crystal_symmetry_operations(:,:,:)
  
  ! Atom sites 
  real(kind=dp), allocatable,     public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp), allocatable,     public, save :: atoms_pos_cart(:,:,:)
  integer, allocatable,           public, save :: atoms_species_num(:)  
  character(len=maxlen), allocatable,  public, save :: atoms_label(:)
  character(len=2), allocatable,  public, save :: atoms_symbol(:)
  integer,                        public, save :: num_atoms
  integer,                        public, save :: num_species

  !-------------------------------------------------------------------------!
  ! G L O B A L L Y   A V A I L A B L E   F U N C T I O N S
  public :: cell_find_MP_grid
  public :: cell_calc_lattice
  public :: cell_report_parameters
  public :: cell_get_atoms
  public :: cell_get_symmetry
  public :: cell_dist
  !-------------------------------------------------------------------------!

contains

  !=========================================================================!
  subroutine cell_find_MP_grid(kpoints,num_kpts,kpoint_grid_dim)
    !=========================================================================!
    ! Find the MP grid from a set of kpoints                                  !
    !-------------------------------------------------------------------------!
    ! Arguments: kpoints - an array of kpoints                                !
    !            num_kpts - size of the kpoint array                          !
    !-------------------------------------------------------------------------!
    ! Returns: kpint_grid_dim - the number of kpoints in each dimension       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: None                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:  None                                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! Described below                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: None                                              
    !--------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------!
    ! Written by Andrew Morris from the lindos program             11/10/2010 !
    !=========================================================================!
    implicit none
    integer, intent(out) :: kpoint_grid_dim(1:3)
    integer, intent(in)  :: num_kpts
    real(kind=dp),intent(in) :: kpoints(1:3,1:num_kpts)
    real(kind=dp) :: kpoints_TR(1:3,1:num_kpts*2)

    ! We have to take into account time reversal symmetry. 
    kpoints_TR(1:3,1:num_kpts)           = kpoints(1:3,1:num_kpts)
    kpoints_TR(1:3,num_kpts+1:num_kpts*2)=-kpoints(1:3,1:num_kpts)

    call kpoint_density(kpoints_TR(1,:), num_kpts*2,  kpoint_grid_dim(1))
    call kpoint_density(kpoints_TR(2,:), num_kpts*2,  kpoint_grid_dim(2))
    call kpoint_density(kpoints_TR(3,:), num_kpts*2,  kpoint_grid_dim(3))
  end subroutine cell_find_MP_grid
  !=========================================================================!

  !=========================================================================!
  subroutine cell_get_symmetry
    !=========================================================================!
    ! Read in the cell symmetries                                             !
    !-------------------------------------------------------------------------!
    ! Arguments: kpoints - an array of kpoints                                !
    !            num_kpts - size of the kpoint array                          !
    !-------------------------------------------------------------------------!
    ! Returns: kpint_grid_dim - the number of kpoints in each dimension       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: None                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:  None                                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! Described below                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: None                                              
    !--------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------!
    ! JRY, April 2011                                                         !
    !=========================================================================!
    use od_comms, only : on_root, comms_bcast
    use od_io, only : filename_len,io_file_unit,seedname,io_error,stdout
    implicit none
    integer :: ierr,sym_file
    logical :: exists
    character(filename_len)     :: sym_filename

    sym_file=io_file_unit()
    sym_filename=trim(seedname)//".sym"
    if(on_root) inquire(file=sym_filename,exist=exists)
    call comms_bcast(exists,1)

    if(.not. exists) then
       write(stdout,*)
       write(stdout,'(80a)') '!--------------------------------- WARNING ------------------------------------!'
       write(stdout,'(80a)') '!                Symmetry Operations file (.sym) not found                     !'
       write(stdout,'(80a)') '!                       Proceeding without symmetry                            !'
       write(stdout,'(80a)') '!------------------------------------------------------------------------------!'
       write(stdout,*)
       num_crystal_symmetry_operations=0
       return
    end if

    if(on_root) then
       open(unit=sym_file,file=sym_filename,form='unformatted',err=100,status='old')

       read(sym_file)num_crystal_symmetry_operations
       if(num_crystal_symmetry_operations > 0)then
          allocate(crystal_symmetry_operations(3,3,num_crystal_symmetry_operations),stat=ierr)
          if(ierr/=0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_get_symmetry")
          allocate(crystal_symmetry_disps(3,num_crystal_symmetry_operations),stat=ierr)
          if(ierr/=0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_get_symmetry")
          read(sym_file)crystal_symmetry_operations
          read(sym_file)crystal_symmetry_disps
       end if
    endif

    call comms_bcast(num_crystal_symmetry_operations,1)
    if(num_crystal_symmetry_operations>0) then
       if(.not. on_root) then
          allocate(crystal_symmetry_operations(3,3,num_crystal_symmetry_operations),stat=ierr)
          if(ierr/=0) call io_error(" Error : cannot allocate crystal_symmetry_operations in cell_get_symmetry")
          allocate(crystal_symmetry_disps(3,num_crystal_symmetry_operations),stat=ierr)
          if(ierr/=0) call io_error(" Error : cannot allocate crystal_symmetry_disps in cell_get_symmetry")
       endif
       call comms_bcast(crystal_symmetry_operations(1,1,1),9*num_crystal_symmetry_operations)
       call comms_bcast(crystal_symmetry_disps(1,1),3*num_crystal_symmetry_operations)
    end if

       return

100 call io_error('Error: Problem opening sym file in cell_get_symmetry') 


  end subroutine cell_get_symmetry

  !=========================================================================
  subroutine kpoint_density(vector,length,points)
    !=========================================================================
    ! Work out the number of kpoints in a given dimension
    ! Receives a vector of length, runs though the vector looking for the
    ! closest distance between points. Then invert this distance to find the
    ! number of k-points in the given dimension. 
    !-------------------------------------------------------------------------
    ! Arguments: vector (in) : vector of kpoints in a particular dimension                              
    !            length (in) : length of vector                  
    !            points (out): number of kpoints in the directon of vector
    !-------------------------------------------------------------------------
    ! Parent module variables used: None                                      
    !-------------------------------------------------------------------------
    ! Modules used:  None                                                     
    !-------------------------------------------------------------------------
    ! Key Internal Variables: Described below                                                         
    !-------------------------------------------------------------------------
    ! Necessary conditions: None            
    !-------------------------------------------------------------------------
    ! Known Worries: None                                      !
    !-------------------------------------------------------------------------
    ! Written by Andrew Morris from the LinDOS program             11/10/2010 
    !=========================================================================
    implicit none
    integer,       intent(in)  :: length
    real(kind=dp), intent(in)  :: vector(1:length)
    integer,       intent(out) :: points

    real(kind=dp) :: real_points
    real(kind=dp) :: distance

    integer :: i,j
    ! THIS IS THE BOMB-PROOF 2nd VERSION
    ! from LinDOS

    ! Check that there isn't only one-kpoint
    if(length==1) then
       points=1
       return
    endif

    ! Check that all of the k-points don't have the same value
    i=1
    do
       if(.not.vector(i)==vector(i+1)) exit
       if((i+1)==length) then ! We haven't exited yet
          ! therefore all must be the same
          points=1
          return
       endif
       i=i+1
    enddo

    ! Now we have more than one k-point and they're not unique
    distance=huge(distance)
    do i=1,length
       do j=i+1,length
          if(vector(i)>vector(j)) then
             if((vector(i)-vector(j))<distance) distance=(vector(i)-vector(j))
          elseif(vector(j)>vector(i)) then
             if((vector(j)-vector(i))<distance) distance=(vector(j)-vector(i))
          endif
       enddo
    enddo


    real_points=1.0_dp/(distance)
    points=int(real_points+0.5_dp)

  end subroutine kpoint_density

  subroutine cell_get_atoms
    use od_constants, only : bohr2ang
    use od_io,        only : io_file_unit, io_error, seedname, maxlen
    use od_algorithms,only : utility_cart_to_frac, utility_frac_to_cart, utility_lowercase

    implicit none

    real(kind=dp),allocatable     :: atoms_pos_frac_tmp(:,:)
    real(kind=dp),allocatable    :: atoms_pos_cart_tmp(:,:)
    character(len=20) :: keyword
    integer           :: in,in1,in2,ins,ine,loop,i,line_e,line_s,counter,tot_num_lines
    integer           :: loop2,max_sites,ierr,ic,num_lines,line_counter,in_unit
    logical           :: found_e,found_s,frac
    character(len=maxlen) :: dummy
    character(len=maxlen),allocatable :: ctemp(:)
    character(len=maxlen),allocatable :: atoms_label_tmp(:)
    logical           :: lconvert

    character(len=maxlen), allocatable :: in_data(:)


    ! read in the cell file

    ! count the lines
    in_unit=io_file_unit( )
    open (in_unit, file=trim(seedname)//'.cell',form='formatted',status='old',err=101)

    num_lines=0;tot_num_lines=0
    do
       read(in_unit, '(a)', iostat = ierr, err= 200, end =210 ) dummy
       dummy=adjustl(dummy)
       tot_num_lines=tot_num_lines+1
       if( .not.dummy(1:1)=='!'  .and. .not. dummy(1:1)=='#' ) then
          if(len(trim(dummy)) > 0 ) num_lines=num_lines+1
       endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.cell')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.cell')
210 continue
    rewind(in_unit)

    ! now read in for real - ignoring comments
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

    ! let's look for the atoms block (remember everything is lower case)
    keyword='positions'


    found_s=.false.
    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'%block')
       if (in==0 .or. in>1) cycle
       if(index(in_data(loop),'frac')>0) then
          frac=.true.
       elseif(index(in_data(loop),'abs')>0) then
          frac=.false.
       else
          cycle
       endif
       line_s=loop
       if (found_s) then
          call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
       endif
       found_s=.true.
    end do

    if(frac) then
       keyword='positions_frac'
    else
       keyword='positions_abs'
    end if

    found_e=.false.
    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'%endblock')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found %block'//trim(keyword)//' more than once in cell file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found %block'//trim(keyword)//' but no %endblock'//trim(keyword)//' in cell file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: %endblock'//trim(keyword)//' comes before %block'//trim(keyword)//' in input file')
    end if

    ! now we know where the atoms block is

    lconvert=.false.
    dummy=in_data(line_s+1)
    if ( index(dummy,'ang').ne.0 ) then
       lconvert=.false.
       line_s=line_s+1
    elseif ( index(dummy,'bohr').ne.0 ) then
       lconvert=.true.
       line_s=line_s+1
    endif

    num_atoms=line_e-1-(line_s+1)+1
    allocate(atoms_pos_frac_tmp(3,num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_frac_tmp in cell_get_atoms')
    allocate(atoms_pos_cart_tmp(3,num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_cart_tmp in cell_get_atoms')
    allocate(ctemp(num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating ctemp in cell_get_atoms')
    allocate(atoms_label_tmp(num_atoms),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_label_tmp in cell_get_atoms')



    counter=0
    do loop=line_s+1,line_e-1
       dummy=in_data(loop)
       counter=counter+1
       if(frac) then
          read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(atoms_pos_frac_tmp(i,counter),i=1,3)
       else
          read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(atoms_pos_cart_tmp(i,counter),i=1,3)
       end if
    end do

    if (lconvert) then
       atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr2ang
    end if


    if(frac) then
       do loop=1,num_atoms
          call utility_frac_to_cart (atoms_pos_frac_tmp(:,loop),atoms_pos_cart_tmp(:,loop),real_lattice)
       end do
    else
       do loop=1,num_atoms
          call utility_cart_to_frac (atoms_pos_cart_tmp(:,loop),atoms_pos_frac_tmp(:,loop),recip_lattice)
       end do
    end if


    ! Now we sort the data into the proper structures
    num_species=1
    ctemp(1)=atoms_label_tmp(1)
    do loop=2,num_atoms
       do loop2=1,loop-1
          if( trim(atoms_label_tmp(loop))==trim(atoms_label_tmp(loop2) )) exit
          if (loop2==loop-1) then 
             num_species=num_species+1
             ctemp(num_species)=atoms_label_tmp(loop)
          end if
       end do
    end do

    allocate(atoms_species_num(num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_species_num in cell_get_atoms')
    allocate(atoms_label(num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_label in cell_get_atoms')
    allocate(atoms_symbol(num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_symbol in cell_get_atoms')
    atoms_species_num(:)=0

    do loop=1,num_species
       atoms_label(loop)=ctemp(loop)
       do loop2=1,num_atoms
          if( trim(atoms_label(loop))==trim(atoms_label_tmp(loop2) )) then
             atoms_species_num(loop)=atoms_species_num(loop)+1
          end if
       end do
    end do

    max_sites=maxval(atoms_species_num)
    allocate(atoms_pos_frac(3,max_sites,num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_frac in cell_get_atoms')
    allocate(atoms_pos_cart(3,max_sites,num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_cart in cell_get_atoms')

    do loop=1,num_species
       counter=0
       do loop2=1,num_atoms
          if( trim(atoms_label(loop))==trim(atoms_label_tmp(loop2) )) then
             counter=counter+1
             atoms_pos_frac(:,counter,loop)=atoms_pos_frac_tmp(:,loop2)
             atoms_pos_cart(:,counter,loop)=atoms_pos_cart_tmp(:,loop2)
          end if
       end do
    end do

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop=1,num_species    
       atoms_symbol(loop)(1:2)=atoms_label(loop)(1:2)
       ic=ichar(atoms_symbol(loop)(2:2))
       if ((ic.lt.ichar('a')).or.(ic.gt.ichar('z'))) &
            atoms_symbol(loop)(2:2)=' '
    end do

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))




  end subroutine cell_get_atoms







  !=========================================================================!
  subroutine cell_calc_lattice
    !=========================================================================!
    ! Begin with a real lattice. Convert from bohr. Calculate a reciprocal
    ! lattice and the volume of the cell
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables used: real_lattice, recip_lattice, cell_volume                                      
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                     
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by Andrew Morris from the LinDOS program             11/10/2010 
    !=========================================================================
    use od_constants, only : pi,bohr2ang
    implicit none

    ! THESE ARE IN BOHR, DON'T GET TRIPPED UP AGAIN!
    real_lattice=real_lattice*bohr2ang

    recip_lattice(1,1)=real_lattice(2,2)*real_lattice(3,3)- &
         real_lattice(3,2)*real_lattice(2,3)
    recip_lattice(2,1)=real_lattice(2,3)*real_lattice(3,1)- &
         real_lattice(3,3)*real_lattice(2,1)
    recip_lattice(3,1)=real_lattice(2,1)*real_lattice(3,2)- &
         real_lattice(3,1)*real_lattice(2,2)
    recip_lattice(1,2)=real_lattice(3,2)*real_lattice(1,3)- &
         real_lattice(1,2)*real_lattice(3,3)
    recip_lattice(2,2)=real_lattice(3,3)*real_lattice(1,1)- &
         real_lattice(1,3)*real_lattice(3,1)
    recip_lattice(3,2)=real_lattice(3,1)*real_lattice(1,2)- &
         real_lattice(1,1)*real_lattice(3,2)
    recip_lattice(1,3)=real_lattice(1,2)*real_lattice(2,3)- &
         real_lattice(2,2)*real_lattice(1,3)
    recip_lattice(2,3)=real_lattice(1,3)*real_lattice(2,1)- &
         real_lattice(2,3)*real_lattice(1,1)
    recip_lattice(3,3)=real_lattice(1,1)*real_lattice(2,2)- &
         real_lattice(2,1)*real_lattice(1,2)

    ! * Calculate cell volume
    cell_volume=real_lattice(1,1)*recip_lattice(1,1) + &
         real_lattice(2,1)*recip_lattice(1,2) + &
         real_lattice(3,1)*recip_lattice(1,3)


    if(cell_volume<0.0_dp) then ! Left handed set
       cell_volume=-cell_volume
    endif

    ! Scale reciprocal lattice by 2*pi/volume
    recip_lattice(:,:)=recip_lattice(:,:)*pi*2.0_dp/cell_volume

  end subroutine cell_calc_lattice

  !=========================================================================!
  subroutine cell_report_parameters
    !=========================================================================
    ! Begin with a real lattice. Convert from bohr. Calculate a reciprocal
    ! lattice and the volume of the cell
    !-------------------------------------------------------------------------
    ! Arguments: None
    !-------------------------------------------------------------------------
    ! Parent module variables used: real_lattice, recip_lattice, cell_volume                                      
    !-------------------------------------------------------------------------
    ! Modules used:  See below                                                 
    !-------------------------------------------------------------------------
    ! Key Internal Variables: None                                                      
    !-------------------------------------------------------------------------
    ! Necessary conditions: None     
    !-------------------------------------------------------------------------
    ! Known Worries: None                                              
    !-------------------------------------------------------------------------
    ! Written by J R Yates, modified A J Morris                     Dec 2010
    !=========================================================================
    use od_io, only : stdout

    implicit none

    integer :: i

    write(stdout,'(30x,a21)') 'Lattice Vectors (Ang)' 
    write(stdout,101) 'a_1',(real_lattice(1,I), i=1,3)
    write(stdout,101) 'a_2',(real_lattice(2,I), i=1,3)
    write(stdout,101) 'a_3',(real_lattice(3,I), i=1,3)

    write(stdout,*)   

    write(stdout,'(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    write(stdout,101) 'b_1',(recip_lattice(1,I), i=1,3)
    write(stdout,101) 'b_2',(recip_lattice(2,I), i=1,3)
    write(stdout,101) 'b_3',(recip_lattice(3,I), i=1,3)

    write(stdout,*) 
    write(stdout,'(19x,a17,3x,f11.5)',advance='no') &
         'Unit Cell Volume:',cell_volume

    write(stdout,'(2x,a7)') '(Ang^3)'
    write(stdout,*) 

    return
101 format(20x,a3,2x,3F11.6)

  end subroutine cell_report_parameters

  subroutine cell_dist
    use od_comms, only : comms_bcast,on_root
    use od_io, only    : io_error
    implicit none
    
    integer :: max_sites,ierr

    call comms_bcast(real_lattice(1,1),9)
    call comms_bcast(recip_lattice(1,1),9)
    call comms_bcast(cell_volume,1)

    !call comms_bcast(kpoint_r(:,:)
    !call comms_bcast(kpoint_r_cart(:,:)
    !call comms_bcast(kpoint_weight(:)


    call comms_bcast(nkpoints,1) 
    call comms_bcast(kpoint_grid_dim(1),3)


    !-------------------------------------------------------------------------!

    call comms_bcast(num_atoms,1)
    call comms_bcast(num_species,1)
    if (num_atoms>0) then
       if(.not. on_root) then
          allocate(atoms_species_num(num_species),stat=ierr)
          if (ierr/=0) call io_error('Error allocating atoms_species_num in cell_dist')
       end if
       call comms_bcast(atoms_species_num(1),num_species)
       max_sites=maxval(atoms_species_num)
       if(.not. on_root) then
          allocate(atoms_pos_frac(3,max_sites,num_species),stat=ierr)
          if (ierr/=0) call io_error('Error allocating atoms_pos_frac in cell_dist')
          allocate(atoms_pos_cart(3,max_sites,num_species),stat=ierr)
          if (ierr/=0) call io_error('Error allocating atoms_pos_cart in cell_dist')
          allocate(atoms_label(num_species),stat=ierr)
          if (ierr/=0) call io_error('Error allocating atoms_label in cell_dist')
          allocate(atoms_symbol(num_species),stat=ierr)
          if (ierr/=0) call io_error('Error allocating atoms_symbol in cell_dist')
       end if
       call comms_bcast(atoms_pos_frac(1,1,1),3*num_species*max_sites)
       call comms_bcast(atoms_pos_cart(1,1,1),3*num_species*max_sites)
       call comms_bcast(atoms_label(1),len(atoms_label(1))*num_species)
       call comms_bcast(atoms_symbol(1),len(atoms_symbol(1))*num_species)
    endif

  end subroutine cell_dist

endmodule od_cell
