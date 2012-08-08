!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! This module contains GPL routines from Wannier90           !
! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          !
!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `COPYING' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!
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

module od_parameters

  use od_constants, only : dp
  use od_io,        only : stdout,maxlen
  use od_cell, only : atoms_label, num_atoms, num_species, atoms_symbol, &
       & atoms_species_num, atoms_pos_frac, atoms_pos_cart, num_crystal_symmetry_operations

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
  logical, public, save :: compare_dos
  logical, public, save :: pdos
  logical, public, save :: jdos
  logical, public, save :: compare_jdos
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
  real(kind=dp),     public, save :: efermi_user     ! If the user has set efermi in the odi file
  character(len=20), public, save :: efermi_choice   ! Where do we want to get the fermi energy from
  logical,           public, save :: finite_bin_correction
  logical,           public, save :: hybrid_linear
  real(kind=dp),     public, save :: hybrid_linear_grad_tol
  logical,           public, save :: numerical_intdos
  logical,           public, save :: compute_band_gap
 
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
  integer, allocatable, public,save :: exclude_bands(:)
  
  integer,           public, save :: num_exclude_bands    ! this is set by param_write

  ! Optics parameters
  character(len=20), public, save :: optics_geom
  real(kind=dp),     public, save :: optics_qdir(3)
  logical,           public, save :: optics_intraband 
  real(kind=dp),     public, save :: optics_drude_broadening 
  real(kind=dp),     public, save :: optics_lossfn_gaussian
  logical,           public, save :: optics_lossfn_broadening ! the is set by param_write


  ! Core parameters 
  character(len=20), public, save :: core_geom
  real(kind=dp),     public, save :: core_qdir(3) 
  logical,           public, save :: core_LAI_broadening
  character(len=20), public, save :: core_type
  real(kind=dp),     public, save :: LAI_gaussian_width
  logical,           public, save :: LAI_gaussian 
  real(kind=dp),     public, save :: LAI_lorentzian_width
  real(kind=dp),     public, save :: LAI_lorentzian_scale 
  real(kind=dp),     public, save :: LAI_lorentzian_offset
  logical,           public, save :: LAI_lorentzian 

  real(kind=dp),     public, save :: lenconfac


  private

  integer :: num_lines
  character(len=maxlen), allocatable :: in_data(:)



  public :: param_read
  public :: param_write
  public :: param_write_header
  public :: param_write_atomic_coord
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
    use od_io,        only : io_error,seedname,stderr
    use od_cell,      only : cell_get_atoms,cell_read_cell

    implicit none

    !local variables
    integer :: i_temp,loop,ierr
    logical :: found
    character(len=20), allocatable :: task_string(:)
    character(len=20) :: c_string



    call param_in_file

    iprint          =  1             ! Verbosity
    call param_get_keyword('iprint',found,i_value=iprint)

    legacy_file_format=.false. 
    call param_get_keyword('legacy_file_format',found,l_value=legacy_file_format)

    energy_unit     =  'ev'          !
    call param_get_keyword('energy_unit',found,c_value=energy_unit)
    if(index(energy_unit,'ev')==0 .and. index(energy_unit,'ry')==0 .and. index(energy_unit,'ha')==0) &
         call io_error('Error: value of energy_unit not recognised in param_read')

    dos=.false.; pdos=.false.; jdos=.false.; optics=.false.; core=.false.; compare_dos=.false.;compare_jdos=.false.
    call param_get_vector_length('task',found,i_temp)
    if(found .and. i_temp>0) then
       allocate(task_string(i_temp),stat=ierr)
       if(ierr/=0) call io_error('Error: param_read - allocation failed for task_string')
       call param_get_keyword_vector('task',found,i_temp,c_value=task_string)
       do loop=1,i_temp
          if(index(task_string(loop),'optics')>0) then
             optics=.true.
          elseif(index(task_string(loop),'core')>0) then
             core=.true.
          elseif(index(task_string(loop),'compare_jdos')>0) then
             jdos=.true.; compare_jdos=.true.
          elseif(index(task_string(loop),'jdos')>0) then
             jdos=.true.
          elseif(index(task_string(loop),'pdos')>0) then
             pdos=.true.
          elseif(index(task_string(loop),'compare_dos')>0) then
             dos=.true.; compare_dos=.true.
          elseif(index(task_string(loop),'dos')>0) then
             dos=.true.
          elseif(index(task_string(loop),'none')>0) then
             dos=.false.; pdos=.false.; jdos=.false.; optics=.false.; core=.false.
          elseif(index(task_string(loop),'all')>0) then
             dos=.true.; pdos=.true.; jdos=.true.; optics=.true.; core=.true.
          else
             call io_error('Error: value of task unrecognised in param_read')
          endif
       end do
       deallocate(task_string,stat=ierr)
       if(ierr/=0) call io_error('Error: param_read - deallocation failed for task_string')
    end if
    if( (compare_dos.or.compare_jdos) .and. (pdos.or.core.or.optics)) &
         call io_error('Error: compare_dos/compare_jdos are not comptable with pdos, core or optics tasks') 

    i_temp=0
    fixed=.false.; adaptive=.false.; linear=.false.; quad=.false. 
    call param_get_keyword('broadening',found,c_value=c_string)
    if (found) then
       if(index(c_string,'fixed')>0) then
          fixed=.true.
       elseif(index(c_string,'adaptive')>0) then
          adaptive=.true.
       elseif(index(c_string,'linear')>0) then
          linear=.true.
       elseif(index(c_string,'quad')>0) then
          quad=.true.
          !          fixed=.true.;adaptive=.true.;linear=.true. 
       else
          call io_error('Error: value of broadening unrecognised in param_read')
       endif
    end if
    if(compare_dos.or.compare_jdos) then
       fixed=.true.;adaptive=.true.;linear=.true.
    end if

    if(.not.(fixed.or.adaptive.or.linear.or.quad)) then ! Pick a default
       adaptive=.true.
    endif

    length_unit     =  'ang'         !
    lenconfac=1.0_dp
    call param_get_keyword('length_unit',found,c_value=length_unit)
    if (length_unit.ne.'ang' .and. length_unit.ne.'bohr') &
         call io_error('Error: value of length_unit not recognised in param_read')
    if (length_unit.eq.'bohr') lenconfac=1.0_dp/bohr2ang

    adaptive_smearing           = 0.4_dp ! LinDOS default
    call param_get_keyword('adaptive_smearing',found,r_value=adaptive_smearing)

    fixed_smearing             = 0.3_dp ! LinDOS default
    call param_get_keyword('fixed_smearing',found,r_value=fixed_smearing)

    efermi_user        = -990.0_dp
    efermi_choice="optados"
    call param_get_efermi('efermi',found,efermi_choice,efermi_user)

    ! Force all Gaussians to be greater than the width of a bin. When using numerical_indos
    ! this is critical for counting all of the Gaussian DOS peaks. 
    ! When using semi-analytic integration it is desirable to show up very sharp peaks in the 
    ! DOS. However, the intDOS will not be affected.
    finite_bin_correction = .true.
    call param_get_keyword('finite_bin_correction',found,l_value=finite_bin_correction)

    ! Perform fixed and adaptive smearing summing the contribution of each Gaussian
    ! instead of the new and better way of taking the erf. Left in for comparison to LinDOS
    numerical_intdos= .false.
    call param_get_keyword('numerical_intdos',found,l_value=numerical_intdos)

    ! Whenever very flat features are found when performing linear broadening, revert to adaptive.
    ! The tolerance is the gradient of the band at the kpoint.
    ! N.B. Finite_bin_correction may also be set, to further improve the spectra.
    hybrid_linear= .false.
    call param_get_keyword('hybrid_linear',found,l_value=hybrid_linear)
    hybrid_linear_grad_tol= 0.01_dp ! Seems about right for getting semi-core states correctly integrated. 
    call param_get_keyword('hybrid_linear_grad_tol',found,r_value=hybrid_linear_grad_tol)
    
    compute_band_energy    = .true.
    call param_get_keyword('compute_band_energy',found,l_value=compute_band_energy)

    set_efermi_zero = .false.
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

    num_exclude_bands=0
    call param_get_range_vector('exclude_bands',found,num_exclude_bands,lcount=.true.)
    if(found) then
       if(num_exclude_bands<1) call io_error('Error: problem reading exclude_bands')
       allocate(exclude_bands(num_exclude_bands),stat=ierr)
       if (ierr/=0) call io_error('Error allocating exclude_bands in param_read')
       call param_get_range_vector('exclude_bands',found,num_exclude_bands,.false.,exclude_bands)
       if (any(exclude_bands<1)  ) &
            call io_error('Error: exclude_bands must contain positive numbers')
    end if

    compute_band_gap        = .false.
    call param_get_keyword('compute_band_gap',found,l_value=compute_band_gap)

    devel_flag=' '
    call param_get_keyword('devel_flag',found,c_value=devel_flag)

    output_format           = 'xmgrace'
    call param_get_keyword('output_format',found,c_value=output_format)
    if(index(output_format,'xmgrace')==0 .and. index(output_format,'gnuplot')==0) &
         call io_error('Error: value of output_format not recognised in param_read')
    

    optics_geom = 'polycrys'
    call param_get_keyword('optics_geom',found,c_value=optics_geom)
    if(index(optics_geom,'polar')==0 .and. index(optics_geom,'polycrys')==0 .and. index(optics_geom,'tensor')==0 ) &
         call io_error('Error: value of optics_geom not recognised in param_read')

    scissor_op          = 0.0_dp !! change
    call param_get_keyword('scissor_op',found,r_value=scissor_op)

    optics_qdir = 0.0_dp
    call  param_get_keyword_vector('optics_qdir',found,3,r_value=optics_qdir)
    if(index(optics_geom,'polar')>0 .and. .not. found) &
         call io_error('Error: polarised or unpolarised optics geometry requested but optics_qdir is not set')
    if( (index(optics_geom,'polycrys')>0 .or. index(optics_geom,'tensor')>0) .and. found) &
         call io_error('Error: polycrystalline optics geometry or full dielectric tensor requested but optics_qdir is set')

    optics_intraband    = .false.
    call param_get_keyword('optics_intraband',found,l_value=optics_intraband)

    optics_drude_broadening          = 1.0e14_dp
    call param_get_keyword('optics_drude_broadening',found,r_value=optics_drude_broadening)
    
    optics_lossfn_broadening=.false.
    optics_lossfn_gaussian=0.0_dp
    call param_get_keyword('optics_lossfn_broadening',optics_lossfn_broadening,r_value=optics_lossfn_gaussian)
    if (optics_lossfn_gaussian<0.0_dp) call io_error('Error: optics_lossfn_broadening must be positive')
    if(abs(optics_lossfn_gaussian)<1.0e-6_dp)  optics_lossfn_broadening=.false. ! trap too small values


    core_geom = 'polycrys'
    call param_get_keyword('core_geom',found,c_value=core_geom)
    if ( (index(core_geom,'polycrys')==0) .and. (index(core_geom,'polar')==0)) &
         call io_error('Error: value of core_geom not recognised in param_read')

    core_type = 'absorption'
    call param_get_keyword('core_type',found,c_value=core_type)
    if (core_type.ne.'absorption' .and. core_type.ne.'emission'.and. core_type.ne.'all') &
         call io_error('Error: value of core_type not recognised in param_read')

    core_qdir = 0.0_dp
    call  param_get_keyword_vector('core_qdir',found,3,r_value=core_qdir)
    if(index(core_geom,'polar')>0 .and. .not. found) &
         call io_error('Error: polarised core geometry requested but core_qdir is not set')
    if(index(core_geom,'polycrys')>0 .and. found) &
         call io_error('Error: polycrystalline core geometry requested but core_qdir is set')

    core_LAI_broadening   = .false.
    call param_get_keyword('core_lai_broadening',found,l_value=core_LAI_broadening)

    LAI_gaussian_width    = 0.0_dp
    LAI_gaussian = .false. 
    call param_get_keyword('lai_gaussian_width',found,r_value=LAI_gaussian_width)
    if (LAI_gaussian_width.gt.1E-14) LAI_gaussian=.true.
    if (LAI_gaussian_width.lt.0.0_dp) call io_error('Error: LAI_gaussian_width must be positive')

    LAI_lorentzian_width    = 0.0_dp
    LAI_lorentzian = .false. 
    call param_get_keyword('lai_lorentzian_width',found,r_value=LAI_lorentzian_width)
    if (LAI_lorentzian_width.gt.1E-14) LAI_lorentzian=.true.
    if (LAI_lorentzian_width.lt.0.0_dp) call io_error('Error: LAI_lorentzian_width must be positive')

    LAI_lorentzian_scale    = 0.1_dp
    call param_get_keyword('lai_lorentzian_scale',found,r_value=LAI_lorentzian_scale)
!    if (LAI_lorentzian_scale.gt.1E-14) LAI_lorentzian=.true. 
    if (LAI_lorentzian_scale.lt.0.0_dp) call io_error('Error: LAI_lorentzian_scale must be positive')

    LAI_lorentzian_offset    = 0.0_dp
    call param_get_keyword('lai_lorentzian_offset',found,r_value=LAI_lorentzian_offset)
    if (LAI_lorentzian_offset.lt.0.0_dp) call io_error('Error: LAI_lorentzian_offset must be positive')

    num_atoms=0
    num_species=0
    num_crystal_symmetry_operations=0
    if(pdos.or.core.or.optics) then
       ! try to read in the atoms from the cell file.
       ! We don't need them otherwise, so let's not bother
       if(index(devel_flag,'old_filename')>0) then
          call cell_get_atoms
       else
          call cell_read_cell
       end if
    end if


    ! check to see that there are no unrecognised keywords

    if ( any(len_trim(in_data(:))>0 )) then
       write(stderr,'(1x,a)') 'The following section of file '//trim(seedname)//'.odi contained unrecognised keywords'
       write(stderr,*) 
       do loop=1,num_lines
          if (len_trim(in_data(loop))>0) then
             write(stderr,'(1x,a)') trim(in_data(loop))
          end if
       end do
       write(stderr,*) 
       call io_error('Unrecognised keyword(s) in input file')
    end if

    call param_uppercase()

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

  end subroutine param_read


  !===================================================================
  subroutine param_uppercase
    !===================================================================
    !                                                                  !
    ! Convert a few things to uppercase to look nice in the output     !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: nsp,ic

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
    write(stdout,'(a78)') " |                       and Jonathan R. Yates                               | "
    write(stdout,'(a78)') " |                                                                           | "
    write(stdout,'(a78)') " |                       Copyright (c) 2010-2012                             | "
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

  subroutine param_write_atomic_coord
    !==================================================================!
    !                                                                  !
    ! write atomic coodes to stdout                                       !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: nat,nsp, atom_counter

    ! System

    if(num_atoms>0) then
       write(stdout,*) ' '
       ! IT DOESN'T SEEM HELPFUL TO WRITE OUT TO MULTIPLY THE INITAL ATOMIC POSITIONS WITH THE FINAL LATTICE...
       if(iprint>2) then
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
          if (lenconfac.eq.1.0_dp) then
             write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
          else
             write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
          endif
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
          do nsp=1,num_species
             do nat=1,atoms_species_num(nsp)
                write(stdout,'(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp),nat,&
                     atoms_pos_frac(:,nat,nsp),'|',atoms_pos_cart(:,nat,nsp)*lenconfac,'|'
             end do
          end do
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
          write(stdout,'(1x,a)') '  WARNING: These are the CASTEP input coordinates not the output -- here to   '
          write(stdout,'(1x,a)') '            aid advanced debugging only.'
       else
          atom_counter=1
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
          write(stdout,'(1x,a)') '|             Species                  Sites                  Total Atoms    |'
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
          do nsp=1,num_species
             write(stdout,'(1x,a1,16x,a2,16x,i4,a3,i4,16x,i4,11x,a)') '|', atoms_symbol(nsp), & 
                  & atom_counter, "to", atom_counter+atoms_species_num(nsp)-1,  atoms_species_num(nsp), "|"
             atom_counter=atom_counter+atoms_species_num(nsp)
          enddo
          write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       endif
       else
       if(iprint>1)  write(stdout,'(25x,a)') 'No atom positions read'
    end if
    write(stdout,*) ' '
  end subroutine param_write_atomic_coord

  !===================================================================
  subroutine param_write
    !==================================================================!
    !                                                                  !
    ! write parameters to stdout                                       !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: nat,nsp

    ! System

   ! if(num_atoms>0) then
   !    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
   !    if (lenconfac.eq.1.0_dp) then
   !       write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
   !    else
   !       write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
   !    endif
   !    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
   !    do nsp=1,num_species
   !       do nat=1,atoms_species_num(nsp)
   !          write(stdout,'(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp),nat,&
   !               atoms_pos_frac(:,nat,nsp),'|',atoms_pos_cart(:,nat,nsp)*lenconfac,'|'
   !       end do
   !    end do
   !    write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
   ! else
   !    if(iprint>1)  write(stdout,'(25x,a)') 'No atom positions read'
   ! end if
  !  write(stdout,*) ' '
!!$
!!$    write(stdout,*) ' '
!!$    if(iprint>1) then
!!$       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
!!$       if (lenconfac.eq.1.0_dp) then
!!$          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
!!$       else
!!$          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
!!$       endif
!!$       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
!!$       !       do nkp=1,nkpoints
!!$       !          write(stdout,'(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',&
!!$       ! nkp,kpoint_r(:,nkp),'|',kpoint_r_cart(:,nkp)/lenconfac,'|'
!!$       !       end do
!!$       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
!!$       write(stdout,*) ' '
!!$    end if
!!$    ! Main


    !
    write(stdout,'(1x,a78)')    '+------------------------------ JOB CONTROL ---------------------------------+'
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
    write(stdout,'(1x,a78)')    '+-------------------------------- UNITS -------------------------------------+'
    write(stdout,'(1x,a46,2x,a4,25x,a1)') '|  Length Unit                               :',trim(length_unit),'|'  

    if(dos.or.pdos) then
       if(dos_per_volume) then 
          write(stdout,'(1x,a78)') '|  J/P/DOS units                             :  electrons eV^-1 Ang^-3       |' 
       else
          write(stdout,'(1x,a78)') '|  J/P/DOS units                             :  electrons eV^-1              |' 
       endif
    endif


    write(stdout,'(1x,a78)')    '+--------------------------SPECTRAL PARAMETERS ------------------------------+'
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
    if(hybrid_linear) then
       write(stdout,'(1x,a78)') '|  Hybrid Linear Correction                     :  True                         |'
       write(stdout,'(1x,a46,2x,F10.8,19x,a1)') '|  Hybrid Linear Gradient Tolerance             :',hybrid_linear_grad_tol,'|'  
    endif
    if(numerical_intdos) &
         write(stdout,'(1x,a78)') '|  Numerical Integration of P/DOS            :  True                         |'        
    if(dos_per_volume) &
         write(stdout,'(1x,a78)') '|  Present DOS per simulation cell volume    :  True                         |'        
    if(set_efermi_zero) then
         write(stdout,'(1x,a78)') '|  Shift energy scale so fermi_energy=0      :  True                         |'        
      else
         write(stdout,'(1x,a78)') '|  Shift energy scale so fermi_energy=0      :  False                        |'        
      end if
    if(compute_band_energy) &
         write(stdout,'(1x,a78)') '|  Compute the band energy                   :  True                         |'        


    if(optics) then
       write(stdout,'(1x,a78)')    '+-------------------------------- OPTICS ------------------------------------+'
       if(index(optics_geom,'polycrys')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Optics Calculation           :  Polycrystalline              |'        
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
       if(optics_intraband) then
          write(stdout,'(1x,a78)') '|  Include Intraband Contribution            :  True                         |'        
          write(stdout,'(1x,a46,1x,1E10.3,20x,a1)') '|  Drude Broadening                          :',&
	& optics_drude_broadening,'|'        
       else
          write(stdout,'(1x,a78)') '|  Include Intraband Contribution            :  False                        |'        
       endif
    end if
    if(core) then
       write(stdout,'(1x,a78)')    '+--------------------------------- CORE -------------------------------------+'
       if(index(core_geom,'polycrys')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Core Calculation             :  Polycrystalline              |'        
       elseif (index(core_geom,'polar')>0) then
          write(stdout,'(1x,a78)') '|  Geometry for Core Calculation             :  Polarised                    |'        
          write(stdout,'(1x,a47,2x,f6.2,2x,f6.2,2x,f6.2,3x,a4)') '|  Direction of q-vector (un-normalised)     : ' &
               ,core_qdir(1:3),'   |'        
       endif
       if(core_LAI_broadening) then
          write(stdout,'(1x,a78)') '|  Include lifetime and Instrument Broadening:  True                         |'        
          write(stdout,'(1x,a46,1x,1f10.4,20x,a1)') '|  Gaussian Width                            :',LAI_gaussian_width,'|'        
          write(stdout,'(1x,a46,1x,1f10.4,20x,a1)') '|  Lorentzian Width                          :',LAI_lorentzian_width,'|'      
          write(stdout,'(1x,a46,1x,1f10.4,20x,a1)') '|  Lorentzian Scale                          :',LAI_lorentzian_scale,'|'
          write(stdout,'(1x,a46,1x,1f10.4,20x,a1)') '|  Lorentzian Offset                         :',LAI_lorentzian_offset,'|'
       else
          write(stdout,'(1x,a78)') '|  Include lifetime and Instrument Broadening:  False                        |'        
       endif
    end if


    write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
    write(stdout,*) ' '

  end subroutine param_write


  !==================================================================!
  subroutine param_dealloc
    !==================================================================!
    !                                                                  !
    ! release memory from allocated parameters                         !
    !                                                                  !
    !===================================================================  

    implicit none

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


  !===========================================================================!
  subroutine param_get_efermi(keyword,found,c_value,r_value)
    !===========================================================================!
    !                                                                           !
    !             Finds the value of the required keyword.                      !
    !                                                                           !
    !===========================================================================!

    use od_io,        only : io_error

    implicit none

    character(*),      intent(in)    :: keyword
    logical,           intent(out)   :: found
    character(*),      intent(inout) :: c_value
    real(kind=dp),     intent(inout) :: r_value

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
       c_value=dummy
       if(trim(c_value)=='optados' .or. trim(c_value)=='file'.or.trim(c_value)=='insulator') then
          r_value=-999.0_dp ! ie not set
       else
          ! assume it is a number
          read(dummy,*,err=220,end=220) r_value
          c_value='user'
       end if
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))


  end subroutine param_get_efermi



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

    use od_comms, only : comms_bcast,on_root
    use od_io   , only : io_error

    implicit none
    integer :: ierr

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
    call comms_bcast(compare_dos  ,1)
    call comms_bcast(compare_jdos ,1)
    call comms_bcast(fixed,1)
    call comms_bcast(adaptive,1)
    call comms_bcast(linear,1)
    call comms_bcast(quad,1)
    call comms_bcast(dos_nbins,1)
    call comms_bcast(compute_band_energy,1)
    call comms_bcast(compute_band_gap,1)
    call comms_bcast(adaptive_smearing,1)
    call comms_bcast(fixed_smearing,1)
    call comms_bcast(hybrid_linear,1)
    call comms_bcast(hybrid_linear_grad_tol,1)
    call comms_bcast(dos_per_volume,1)
    call comms_bcast(efermi_user,1)
    call comms_bcast(efermi_choice,len(efermi_choice))
    call comms_bcast(finite_bin_correction,1)
    call comms_bcast(numerical_intdos ,1)
    call comms_bcast(jdos_max_energy ,1)
    call comms_bcast(jdos_spacing ,1)
    call comms_bcast(scissor_op,1)
    call comms_bcast(optics_geom,len(optics_geom))
    call comms_bcast(optics_qdir(1),3)
    call comms_bcast(optics_intraband,1)
    call comms_bcast(optics_drude_broadening,1)
    call comms_bcast(core_geom,len(core_geom))
    call comms_bcast(core_type,len(core_type))
    call comms_bcast(core_qdir(1),3)
    call comms_bcast(core_LAI_broadening,1)
    call comms_bcast(LAI_gaussian_width,1)
    call comms_bcast(LAI_gaussian,1)
    call comms_bcast(LAI_lorentzian_width,1)
    call comms_bcast(LAI_lorentzian_scale,1)
    call comms_bcast(LAI_lorentzian_offset,1)
    call comms_bcast(LAI_lorentzian,1)
    call comms_bcast(dos_per_volume,1)
    call comms_bcast(dos_min_energy,1)
    call comms_bcast(dos_max_energy,1)
    call comms_bcast(dos_spacing,1)
    call comms_bcast(legacy_file_format,1)
    call comms_bcast(pdos_string,len(pdos_string))
    call comms_bcast(set_efermi_zero,1)
    !
    call comms_bcast(num_exclude_bands,1)
    if(num_exclude_bands>1) then
       if(.not.on_root) then
          allocate(exclude_bands(num_exclude_bands),stat=ierr)
          if (ierr/=0) call io_error('Error allocating exclude_bands in param_read')
       endif
       call comms_bcast(exclude_bands(1),num_exclude_bands)
    end if


  end subroutine param_dist

end module od_parameters
