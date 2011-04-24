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

  use od_constants, only : dp
  use od_io,        only : stdout,maxlen
  use od_cell, only : real_lattice, recip_lattice, cell_volume, kpoint_grid_dim, nkpoints, &
       & kpoint_r, kpoint_r_cart, atoms_label, num_atoms, num_species, atoms_symbol, &
       & atoms_species_num, atoms_pos_frac, atoms_pos_cart

  implicit none

  !Generic Parameters
  character(len=20), public, save :: output_format
  character(len=100), public, save :: devel_flag
  integer,           public, save :: iprint
  character(len=20), public, save :: energy_unit
  character(len=20), public, save :: length_unit ! not exposed - but useful for BZ plots?
  logical,           public, save :: legacy_file_format


  !Task parameters
  logical, public, save :: dos
  logical, public, save :: pdos
  logical, public, save :: jdos
  logical, public, save :: optics
  logical, public, save :: core

  !Broadening parameters
  logical, public, save :: fixed
  logical, public, save :: adaptive
  logical, public, save :: linear
  logical, public, save :: quad

  ! Belonging to the dos module
  logical, public, save :: compute_band_energy
  real(kind=dp),     public, save :: adaptive_smearing
  real(kind=dp),     public, save :: fixed_smearing
  logical,           public, save :: dos_per_volume
  real(kind=dp),     public, save :: fermi_energy
  logical,           public, save :: compute_efermi
  logical,           public, save :: finite_bin_correction
  logical,           public, save :: numerical_intdos 

  logical,           public, save :: set_efermi_zero  ! Set fermi level to zero in dos plot default true
  real(kind=dp),     public, save :: dos_min_energy
  real(kind=dp),     public, save :: dos_max_energy
  real(kind=dp),     public, save :: dos_spacing
  integer,           public, save :: dos_nbins

  ! pdos
  character(len=maxlen), public, save :: pdos_string


  ! Belonging to the jdos module
  real(kind=dp),     public, save :: jdos_max_energy 
  real(kind=dp),     public, save :: jdos_spacing
  real(kind=dp),     public, save :: scissor_op

  ! Optics parameters
  character(len=20), public, save :: optics_geom
  real(kind=dp),     public, save :: optics_qdir(3)

  real(kind=dp),     public, save :: lenconfac

  private

  integer :: num_lines
  character(len=maxlen), allocatable :: in_data(:)



  public :: param_read
  public :: param_write
  public :: param_write_header
  public :: param_dealloc
  public :: param_dist


contains

  !==================================================================!
  subroutine param_read ( )
    !==================================================================!
    !                                                                  !
    ! Read parameters and calculate derived values                     !
    !                                                                  !
    !===================================================================  

    use od_constants, only : bohr2ang
    use od_io,        only : io_error,io_file_unit,seedname
    use od_cell,      only : cell_get_atoms
!    use od_electronic,only : elec_pdos_read_orbitals
    implicit none

    !local variables
    real(kind=dp)  :: real_lattice_tmp(3,3)
    integer :: nkp,i,j,n,k,itmp,i_temp,i_temp2,eig_unit,loop,ierr,iv_temp(3)
    logical :: found,found2,eig_found,lunits,chk_found
    character(len=6) :: spin_str
    real(kind=dp) :: cosa(3),rv_temp(3)
    character(len=10), allocatable :: task_string(:)

    

    call param_in_file

    iprint          =  1             ! Verbosity
    call param_get_keyword('iprint',found,i_value=iprint)

    legacy_file_format=.true. 
    call param_get_keyword('legacy_file_format',found,l_value=legacy_file_format)

    energy_unit     =  'ev'          !
    call param_get_keyword('energy_unit',found,c_value=energy_unit)

    dos=.false.; pdos=.false.; jdos=.false.; optics=.false.; core=.false.
    call param_get_vector_length('task',found,i_temp)
    if(found .and. i_temp>0) then
       allocate(task_string(i_temp))
       call param_get_keyword_vector('task',found,i_temp,c_value=task_string)
       do loop=1,i_temp
          if(index(task_string(loop),'optics')>0) then
             optics=.true.
          elseif(index(task_string(loop),'core')>0) then
             core=.true.
          elseif(index(task_string(loop),'jdos')>0) then
             jdos=.true.
          elseif(index(task_string(loop),'pdos')>0) then
             pdos=.true.
          elseif(index(task_string(loop),'dos')>0) then
             dos=.true.
          elseif(index(task_string(loop),'none')>0) then
             dos=.false.; pdos=.false.; jdos=.false.; optics=.false.; core=.false.
          elseif(index(task_string(loop),'all')>0) then
             dos=.true.; pdos=.false.; jdos=.false.; optics=.false.; core=.false.
          else
             call io_error('Error: value of task unrecognised in param_read')
          endif
       end do
       deallocate(task_string)
    end if

    num_atoms=0
    num_species=0
    if(pdos) then
       ! try to read in the atoms from the cell file.
       ! We don't need them otherwise, so let's not bother
       call cell_get_atoms
    end if

    i_temp=0
    fixed=.false.; adaptive=.false.; linear=.false.; quad=.false. 
    call param_get_vector_length('broadening',found,i_temp)
    if(found .and. i_temp>0) then
       allocate(task_string(i_temp))
       call param_get_keyword_vector('broadening',found,i_temp,c_value=task_string)
       do loop=1,i_temp
          if(index(task_string(loop),'fixed')>0) then
             fixed=.true.
          elseif(index(task_string(loop),'adaptive')>0) then
             adaptive=.true.
          elseif(index(task_string(loop),'linear')>0) then
             linear=.true.
          elseif(index(task_string(loop),'quad')>0) then
             quad=.true.
          elseif(index(task_string(loop),'all')>0) then
             fixed=.true.;adaptive=.true.;linear=.true. 
          else
             call io_error('Error: value of broadening unrecognised in param_read')
          endif
       end do
       deallocate(task_string)
    end if

    if(.not.(fixed.or.adaptive.or.linear.or.quad)) then ! Piak a default
       adaptive=.true.
    endif

    length_unit     =  'ang'         !
    lenconfac=1.0_dp
    call param_get_keyword('length_unit',found,c_value=length_unit)
    if (length_unit.ne.'ang' .and. length_unit.ne.'bohr') &
         call io_error('Error: value of length_unit not recognised in param_read')
    if (length_unit.eq.'bohr') lenconfac=1.0_dp/bohr2ang

    adaptive_smearing           = 1.4_dp ! LinDOS default
    call param_get_keyword('adaptive_smearing',found,r_value=adaptive_smearing)

    fixed_smearing             = 0.3_dp ! LinDOS default
    call param_get_keyword('fixed_smearing',found,r_value=fixed_smearing)

    fermi_energy        = -990.0_dp
    call param_get_keyword('fermi_energy',found,r_value=fermi_energy)

    ! Force all Gaussians to be greater than the width of a bin. When using numerical_indos
    ! this is critical for counting all of the Gaussian DOS peaks. 
    ! When using semi-analytic integration it is desirable to show up very sharp peaks in the 
    ! DOS. However, the intDOS will not be affected.
    finite_bin_correction = .false.
    call param_get_keyword('finite_bin_correction',found,l_value=finite_bin_correction)

    ! Perform fixed and adaptive smearing summing the contribution of each Gaussian
    ! instead of the new and better way of taking the erf. Left in for comparison to LinDOS
    numerical_intdos= .false.
    call param_get_keyword('numerical_intdos',found,l_value=numerical_intdos)

    compute_band_energy    = .true.
    call param_get_keyword('compute_band_energy',found,l_value=compute_band_energy)

    set_efermi_zero = .true. 
    call param_get_keyword('set_efermi_zero',found,l_value=set_efermi_zero)

    dos_per_volume = .false.
    call param_get_keyword('dos_per_volume',found,l_value=dos_per_volume)

    dos_min_energy         = -huge(dos_min_energy)
    call param_get_keyword('dos_min_energy',found,r_value=dos_min_energy)

    dos_max_energy         = huge(dos_max_energy)
    call param_get_keyword('dos_max_energy',found,r_value=dos_max_energy)

    dos_spacing            = -1.0_dp
    call param_get_keyword('dos_spacing',found,r_value=dos_spacing)

    dos_nbins               = -1 ! 10001 LinDOS default
    call param_get_keyword('dos_nbins',found,i_value=dos_nbins)

    pdos_string =''
    call param_get_keyword('pdos',found,c_value=pdos_string)
    if(pdos.and. (len_trim(pdos_string)==0)) call io_error('pdos requested but pdos is not specified')

    jdos_max_energy        = -1.0_dp !! change
    call param_get_keyword('jdos_max_energy',found,r_value=jdos_max_energy)

    jdos_spacing          = 0.01_dp !! change
    call param_get_keyword('jdos_spacing',found,r_value=jdos_spacing)

    compute_efermi        = .false.
    call param_get_keyword('compute_efermi',found,l_value=compute_efermi)

    devel_flag=' '
    call param_get_keyword('devel_flag',found,c_value=devel_flag)

    output_format           = 'xmgrace'
    call param_get_keyword('output_format',found,c_value=output_format)

    optics_geom = 'polycrys'
    call param_get_keyword('optics_geom',found,c_value=optics_geom)

    if(index(optics_geom,'polar')==0 .and. index(optics_geom,'polycrys')==0 .and. index(optics_geom,'tensor')==0 ) &
         call io_error('Error: value of optics_geom not recognised in param_read')

    optics_qdir = 0.0_dp
    call  param_get_keyword_vector('optics_qdir',found,3,r_value=optics_qdir)
    if(index(optics_geom,'polar')>0 .and. .not. found) &
         call io_error('Error: polarised or unpolarised optics geometry requested but optics_qdir is not set')
    if( (index(optics_geom,'polycrys')>0 .or. index(optics_geom,'tensor')>0) .and. found) &
         call io_error('Error: polycrystalline optics geometry or full dielectric tensor requested but optics_qdir is set')


    call param_uppercase()

303 continue

    deallocate(in_data,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating in_data in param_read')

    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !


    if(dos_min_energy.ge.dos_max_energy) then
       call io_error('Error: must have dos_min_energy < dos_max_energy')
    endif

    if((dos_nbins>0).and.(dos_spacing>0.0_dp)) then
       call io_error('Error: only one of dos_nbins and dos_spacing may be set')
    endif

    if((dos_nbins<0).and.(dos_spacing<0.0_dp)) then
       dos_spacing= 0.005 ! Roughly similar to LinDOS 
    endif


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


    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic=ichar(length_unit(1:1))
    if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
         length_unit(1:1) = char(ic+ichar('Z')-ichar('z'))

    return

  end subroutine param_uppercase


  !===================================================================
  subroutine param_write_header
    use od_constants, only :  optados_version
    implicit none
    write(stdout,*)
    write(stdout,'(a78)') " +===========================================================================+"
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |                OOO   PPPP  TTTTT  AA   DDD    OOO    SSS                  | "
    write(stdout,'(a78)') " |               O   O  P   P   T   A  A  D  D  O   O  S                     | "
    write(stdout,'(a78)') " |               O   O  PPPP    T   AAAA  D  D  O   O   SS                   | "
    write(stdout,'(a78)') " |               O   O  P       T   A  A  D  D  O   O     S                  | "
    write(stdout,'(a78)') " |                OOO   P       T   A  A  DDD    OOO   SSS                   | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " +---------------------------------------------------------------------------+ "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a46,a5,a28)') " |                 Welcome to OptaDOS version ", optados_version,"   &
         &                     | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |         Andrew J. Morris, Rebecca Nicholls, Chris J. Pickard              | "
    write(stdout,'(a78)') " |                       and Jonathan Yates                                  | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |                       Copyright (c) 2010                                  | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |  Please cite:                                                             | "
    write(stdout,'(a78)') " |  Andrew J. Morris, Rebecca Nicholls, Chris J. Pickard and Jonathan Yates  | "
    write(stdout,'(a78)') " |    OptaDOS User Manual, Univ. of Oxford and Univ. College London, (2010)  | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |  Additionally when using the linear broadening:                           | "
    write(stdout,'(a78)') " | C.J. Pickard and M.C. Payne, Phys. Rev. B, 59, 7, 4685 (1999)             | "
    write(stdout,'(a78)') " | C.J. Pickard and M.C. Payne, Phys. Rev. B, 62, 7, 4383 (2000)             | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |  Additionally when using the adaptive broadening:                         | "
    write(stdout,'(a78)') " | J.Yates, X.Wang, D.Vanderbilt and I.Souza, Phys. Rev. B, 75, 195121 (2007)| "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |      in all your publications arising from your use of OptaDOS            | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " +===========================================================================+ "
    write(stdout,*)
  end subroutine param_write_header


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
       if(iprint>1)  write(stdout,'(25x,a)') 'No atom positions read'
    end if
    write(stdout,*) ' '

    write(stdout,*) ' '
    if(iprint>1) then
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
       else
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
       endif
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       !       do nkp=1,nkpoints
       !          write(stdout,'(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',nkp,kpoint_r(:,nkp),'|',kpoint_r_cart(:,nkp)/lenconfac,'|'
       !       end do
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       write(stdout,*) ' '
    end if
    ! Main


    !
    write(stdout,'(1x,a78)')    '*-------------------------------- TASK --------------------------------------*'
    !
    if(dos) then
       write(stdout,'(1x,a78)') '|  Output Density of States                  :  True                         |'
    else
       write(stdout,'(1x,a78)') '|  Output Density of States                  :  False                        |'
    endif
    if(pdos) then
       write(stdout,'(1x,a78)') '|  Output Partial Density of States          :  True                         |'
    else
       write(stdout,'(1x,a78)') '|  Output Partial Density of States          :  False                        |'
    endif
    if(jdos) then
       write(stdout,'(1x,a78)') '|  Output Joint Density of States            :  True                         |'
    else
       write(stdout,'(1x,a78)') '|  Output Joint Density of States            :  False                        |'
    endif
    if(optics) then
       write(stdout,'(1x,a78)') '|  Output Optical Response                   :  True                         |'
    else
       write(stdout,'(1x,a78)') '|  Output Optical Response                   :  False                        |'
    endif
    if(core) then
       write(stdout,'(1x,a78)') '|  Output Core-level Spectra                 :  True                         |'
    else
       write(stdout,'(1x,a78)') '|  Output Core-level Spectra                 :  False                        |'
    endif
    write(stdout,'(1x,a78)')    '*-------------------------------- UNITS -------------------------------------*'
    write(stdout,'(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :',trim(length_unit),'|'  

    if(dos.or.pdos) then
       if(dos_per_volume) then 
          write(stdout,'(1x,a78)') '|  J/P/DOS units                             :  electrons eV^-1 Ang^-3       |' 
       else
          write(stdout,'(1x,a78)') '|  J/P/DOS units                             :  electrons eV^-1              |' 
       endif
    endif


    write(stdout,'(1x,a78)')    '*-------------------------------- BROADENING --------------------------------*'
    if(fixed) then
       write(stdout,'(1x,a78)') '|  Fixed Width Smearing                      :  True                         |'
       write(stdout,'(1x,a46,1x,1F10.5,20x,a1)') '|  Smearing Width                            :', fixed_smearing,'|'
    endif
    if(adaptive) then
       write(stdout,'(1x,a78)') '|  Adaptive Width Smearing                   :  True                         |'
       write(stdout,'(1x,a46,1x,1F10.5,20x,a1)') '|  Adaptive Smearing ratio                   :', adaptive_smearing,'|'
    endif
    if(linear) &
         write(stdout,'(1x,a78)') '|  Linear Extrapolation                      :  True                         |'
    if(quad) &
         write(stdout,'(1x,a78)') '|  Quadratic Extrapolation                   :  True                         |'
    if(finite_bin_correction) &
         write(stdout,'(1x,a78)') '|  Finite Bin Correction                     :  True                         |'
    if(numerical_intdos) &
         write(stdout,'(1x,a78)') '|  Numerical Integration of P/DOS            :  True                         |'        


!    write(stdout,'(1x,a78)')    '*----------------------------- Parameters -----------------------------------*'
!    write(stdout,'(1x,a78)')    '|                                                                            |'
    if(optics) then
       write(stdout,'(1x,a78)')    '*-------------------------------- OPTICS ------------------------------------*'
       if(index(optics_geom,'polycrys')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Optics Calculation           :  Polycrytalline               |'        
       elseif (index(optics_geom,'unpolar')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Optics Calculation           :  Unpolarised                  |'        
          write(stdout,'(1x,a47,2x,f6.2,2x,f6.2,2x,f6.2,3x,a4)') '|  Direction of q-vector (un-normalised)     : ' &
               ,optics_qdir(1:3),'   |'        
       elseif (index(optics_geom,'polar')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Optics Calculation           :  Polarised                    |'        
          write(stdout,'(1x,a47,2x,f6.2,2x,f6.2,2x,f6.2,3x,a4)') '|  Direction of q-vector (un-normalised)     : ' &
               ,optics_qdir(1:3),'   |'        
       elseif (index(optics_geom,'tensor')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Optics Calculation           :  Full dielectric tensor       |'        
       end if
    end if
    write(stdout,'(1x,a78)')    '------------------------------------------------------------------------------'


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
    !       ieallocate (  ndimwin, stat=ierr  )
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
    use od_algorithms,only : utility_lowercase

    implicit none

    integer           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy

    in_unit=io_file_unit( )
    open (in_unit, file=trim(seedname)//'.odi',form='formatted',status='old',err=101)

    num_lines=0;tot_num_lines=0
    do
       read(in_unit, '(a)', iostat = ierr, err= 200, end =210 ) dummy
       dummy=adjustl(dummy)
       tot_num_lines=tot_num_lines+1
       if( .not.dummy(1:1)=='!'  .and. .not. dummy(1:1)=='#' ) then
          if(len(trim(dummy)) > 0 ) num_lines=num_lines+1
       endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.odi')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.odi')
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

  subroutine param_get_pdos

    !    call param_get_block_length('projection',found,rows)




  end subroutine param_get_pdos


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

    use od_constants, only : bohr2ang
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
          r_value=r_value*bohr2ang
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


  subroutine param_dist
    !-----------------------------------------------------
    ! Send the parameters from the root node to all others
    !-----------------------------------------------------

    use od_comms, only : comms_bcast

    implicit none

    call comms_bcast(output_format,len(output_format))
    call comms_bcast(devel_flag   ,len(devel_flag))
    call comms_bcast(iprint        ,1)
    call comms_bcast(energy_unit   ,len(energy_unit))
    call comms_bcast(length_unit ,len(length_unit))
    call comms_bcast(dos    ,1)
    call comms_bcast(pdos   ,1)
    call comms_bcast(jdos   ,1)
    call comms_bcast(optics ,1)
    call comms_bcast(core   ,1)
    call comms_bcast(fixed,1)
    call comms_bcast(adaptive,1)
    call comms_bcast(linear,1)
    call comms_bcast(quad,1)
    call comms_bcast(dos_nbins,1)
    call comms_bcast(compute_band_energy,1)
    call comms_bcast(adaptive_smearing,1)
    call comms_bcast(fixed_smearing,1)
    call comms_bcast(dos_per_volume,1)
    call comms_bcast(fermi_energy,1)
    call comms_bcast(compute_efermi,1)
    call comms_bcast(finite_bin_correction,1)
    call comms_bcast(numerical_intdos ,1)
    call comms_bcast(jdos_max_energy ,1)
    call comms_bcast(jdos_spacing ,1)
    call comms_bcast(scissor_op,1)
    call comms_bcast(optics_geom,len(optics_geom))
    call comms_bcast(optics_qdir(1),3)
    call comms_bcast(dos_per_volume,1)
    call comms_bcast(dos_min_energy,1)
    call comms_bcast(dos_max_energy,1)
    call comms_bcast(dos_spacing,1)
    call comms_bcast(legacy_file_format,1)
    call comms_bcast(pdos_string,len(pdos_string))

  end subroutine param_dist

end module od_parameters
