!=============================================================================== 
! MODULE od_jdos  OptaDOS - Joint Density of States
! This is the module that contains all if the JDOS routines. It is used through 
! the global jdos_calculate subroutine
!-------------------------------------------------------------------------------
! NB. It should be possible to pass optioinal arguments to sub programs as
! optional argumnets without checking whether they are there or not. g95 will 
! allow this behaviour. gfotran will not.
!=============================================================================== 
module od_jdos
 use od_algorithms, only : heap_sort,gaussian
 use od_constants,  only : bohr, dp, H2eV
 use od_comms,      only : on_root, num_nodes, my_node_id, root_id,comms_slice,comms_bcast,& 
         &comms_send, comms_recv,comms_reduce
 use od_electronic, only : band_energy,band_gradient,efermi,efermi_castep,num_electrons, &
         &nbands, nspins, spin_polarised, elec_report_parameters, elec_read_band_gradient, &
         &elec_read_band_energy, electrons_per_state
 use od_io,         only : stdout,io_date,io_error,maxlen,seedname,filename_len, &
         &io_file_unit,io_time
 use od_cell,       only : real_lattice, recip_lattice, cell_volume, cell_find_MP_grid, & 
         &kpoint_r, kpoint_weight, nkpoints, kpoint_grid_dim, cell_calc_lattice, &
         &cell_report_parameters
 use od_dos_utils,        only : dos_utils_merge, doslin_sub_cell_corners, doslin
 use od_parameters

 implicit none


!-------------------------------------------------------------------------------
! P U B L I C   V A R I A B L E S
 real(kind=dp), allocatable, public, save :: jdos_adaptive(:,:)
 real(kind=dp), allocatable, public, save :: jdos_fixed(:,:)
 real(kind=dp), allocatable, public, save :: jdos_linear(:,:)

 real(kind=dp), allocatable, public, save :: E(:)
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! P U B L I C   F U N C T I O N S 
 public :: jdos_calculate
!-------------------------------------------------------------------------------

 real(kind=dp), save                   :: delta_bins ! Width of bins
 integer, allocatable, save :: vb_max(:)
!-------------------------------------------------------------------------------

contains

!=============================================================================== 
subroutine jdos_calculate
!=============================================================================== 
! Main routine in dos module, drives the calculation of Density of states for
! both task : dos and also if it is required elsewhere.
!=============================================================================== 
implicit none
integer :: ierr, idos, i, ik, is, ib
real(kind=dp) :: time0, time1

if(on_root) then
  write(stdout,*)
  write(stdout,'(1x,a78)')'+============================================================================+'
  write(stdout,'(1x,a78)')'+================== Joint Density Of States Calculation =====================+'
  write(stdout,'(1x,a78)')'+============================================================================+'
  write(stdout,*)
endif

!-------------------------------------------------------------------------------
! R E A D   B A N D S   F I L E
! The .band file contains a lot of other important information then just the bands.
! We cannot write this out prior to the .bands file being read.

call elec_read_band_energy

if(on_root) call cell_calc_lattice
if(on_root) call cell_report_parameters
if(on_root) call elec_report_parameters
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! R E A D   B A N D   G R A D I E N T S 
! If we're using one of the more accurate roadening schemes we also need to read in the 
! band gradients too
if(quad.or.linear.or.adaptive) call elec_read_band_gradient
!-------------------------------------------------------------------------------


allocate(vb_max(nspins), stat=ierr)
if (ierr/=0) call io_error ("cannot allocate vb_max")
  
  
! For an insulator
vb_max(:)=num_electrons(:)/electrons_per_state


!-------------------------------------------------------------------------------
! C A L C U L A T E   D O S 
! Now everything is set up, we can perform the dos accumulation in parellel
time0=io_time()

call setup_energy_scale
    
if(on_root.and.(iprint>1)) write(stdout,*)

if(fixed)then
     call calculate_jdos(jdos_fixed)
 
    call dos_utils_merge(jdos_fixed) 
  
 endif
if(adaptive)then
    call calculate_jdos(jdos_adaptive)

    call dos_utils_merge(jdos_adaptive)

  endif
if(linear)then
    call calculate_jdos(jdos_linear)
  
    call dos_utils_merge(jdos_linear)
endif


if(quad) then
  call io_error("quadratic broadening not implemented") 
 !if(quad)    call merge_dos(dos_quad)
 !if(quad)    call merge_dos(intdos_quad)
endif  
  
if(.not.on_root) then
  if(allocated(E)) deallocate(E, stat=ierr)
  if (ierr/=0) call io_error ("cannot deallocate  E")
endif

time1=io_time()
if(on_root)  write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate dos  ',time1-time0,' (sec)'
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! F E R M I   E N E R G Y   A N A L Y S I S
!time0=io_time()
!write(stdout,*)
!write(stdout,'(1x,a78)')  '+------------------------ Fermi Energy Analysis -----------------------------+'
!write(stdout,'(1x,a1,a45,f8.4,a3,20x,a1)') "|","Fermi energy from CASTEP :",efermi_castep," eV","|"
!write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
!
!if(compute_efermi) then
! if(fixed) then 
!   write(stdout,'(1x,a78)') "| From Fixed broadening                                                      | "
!   efermi_fixed= calc_efermi_from_intdos(intdos_fixed)
!   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Fixed braodening) : ", Efermi_fixed,"eV","|"
!   write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'

! endif
! if(adaptive) then
!   write(stdout,'(1x,a78)') "| From Adaptive broadening                                                   | "  
!   efermi_adaptive=calc_efermi_from_intdos(intdos_adaptive)
!   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Adaptive braodening) : ", Efermi_adaptive,"eV","|"
! write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'

! endif
! if(linear) then
!   write(stdout,'(1x,a78)') "| From Linear broadening                                                     | " 
!   efermi_linear=calc_efermi_from_intdos(intdos_linear) 
!   write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Linear braodening) : ", Efermi_linear," eV","|"
! write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
!
! endif
!else
! write(stdout,'(1x,a78)')   "| No Fermi energies calculated                                               | "
! ! Use the derived Fermi energy shifts if calculated. It not use the fermi_energy supplied
! ! if not, use the CASTEP Fermi energy
! if(fermi_energy.ne.-990.0_dp)then
!  efermi_fixed=fermi_energy
!  efermi_adaptive=fermi_energy
!  efermi_linear=fermi_energy
! else
!  efermi_fixed=efermi_castep
!  efermi_adaptive=efermi_castep
!  efermi_linear=efermi_castep
! endif
!endif

! NB If you have asked for more than one type of broadening
! If one of your options is linear then all will have the linear efermi
! If you have asked for adaptive, but nor linear, you will have the adaptive efermi
!if(fixed) then
! efermi=efermi_fixed  
!endif

!if(adaptive) then
! efermi=efermi_adaptive
!endif
!
!if(linear)then
! efermi=efermi_linear
!endif
!
! write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy used : ", efermi,"eV","|"
! E(:)=E(:)-efermi
! band_energy(:,:,:) = band_energy(:,:,:) - efermi
!
!write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
!time1=io_time()
!write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Fermi energies ',time1-time0,' (sec)'
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! B A N D   E N E R G Y   A N A L Y S I S
! Now for a bit of crosschecking  band energies
! These should all converge to the same number as the number of bins is increased
!if(compute_band_energy) call compute_band_energies
!-------------------------------------------------------------------------------


if(dos_per_volume) then
 if(fixed) then
    jdos_fixed=jdos_fixed/cell_volume   
 endif
 if(adaptive) then
    jdos_adaptive=jdos_adaptive/cell_volume 
 endif
 if(linear) then
     jdos_linear=jdos_linear/cell_volume     
 endif

! if(quad) then
!    dos_quad=dos_quad/cell_volume         
!    intdos_quad=intdos_quad/cell_volume
! endif   
endif

!-------------------------------------------------------------------------------
! W R I T E   O U T   D O S  
if(jdos) then ! We have to write stuff out
time0=io_time()
! Otherwise we have written to wdos and dos, so they can be called 
! by whatever.
  if(on_root) then
    if(fixed)    call write_dos(E, jdos_fixed,  "fixed")
    if(adaptive) call write_dos(E, jdos_adaptive, "adaptive")
    if(linear)   call write_dos(E, jdos_linear,  "linear")
    !if(quad)    call write_dos(E, dos_quad, intdos_quad, "quad")
   endif
time1=io_time()
write(stdout,'(1x,a40,f11.3,a)') 'Time to write dos to disk ',time1-time0,' (sec)'
else
 write(stdout,'(1x,a40)') 'Skipping writing out DOS to file'
endif
!-------------------------------------------------------------------------------


write(stdout,'(1x,a78)')    '+============================================================================+'
write(stdout,'(1x,a78)')    '+============== Joint Density Of States Calculation End =====================+'
write(stdout,'(1x,a78)')    '+============================================================================+'
write(stdout,*)

end subroutine jdos_calculate


!=============================================================================== 
subroutine write_dos(E,dos,dos_name)
!=============================================================================== 
! This routine receives an energy scale, a density of states and a file name
! and writes out the DOS to disk
!=============================================================================== 
  implicit none
  real(dp), intent(in) :: E(nbins)
  real(dp), intent(in) :: dos(nbins,nspins)
  character(len=*), intent(in) :: dos_name
  integer :: i, dos_file, ierr
  character(len=11) :: cdate
  character(len=9) :: ctime
  character(len=20) :: dos_units, intdos_units
  

  dos_file=io_file_unit()
  open(unit=dos_file,file=trim(seedname)//'.j'//trim(dos_name)//'.dat',iostat=ierr)
  if(ierr.ne.0) call io_error(" ERROR: Cannot open output file in dos: write_dos")

  dos_units="(electrons per eV)" ; intdos_units="(electrons)"
  if(dos_per_volume) then
    dos_units="(electrons per eV/A^3)" 
    intdos_units="(electrons per A^3)"
  endif

  write(dos_file, *) "##############################################################################"
  write(dos_file,*) "#"
  write(dos_file, *) "#                  O p t a D O S   o u t p u t   f i l e "  
  write(dos_file, '(1x,a1)') "#"
  write(dos_file,*) "#    Denisty of States using ", trim(dos_name), " broadening"
  call io_date(cdate,ctime)
  write(dos_file,*)  '#  Generated on ',cdate,' at ',ctime
  write(dos_file,*) "# Column        Data"
  write(dos_file,*) "#    1        Energy (eV)"
  if(nspins>1) then
    write(dos_file,*) "#    1        Up-spin DOS ", trim(dos_units)
    write(dos_file,*) "#    2        D50own-spin DOS ", trim(dos_units)
    write(dos_file,*) "#    3        Up-spin Integrated DOS ", trim(intdos_units)
    write(dos_file,*) "#    4        Down-spin Integrated DOS ", trim(intdos_units)
  else
    write(dos_file,*) "#    2        DOS ", trim(dos_units)
    write(dos_file,*) "#    3        Integrated DOS ", trim(intdos_units)
  endif
  write(dos_file, '(1x,a1)') "#"
  write(dos_file, '(1x,a78)') "##############################################################################"


  if(nspins>1) then
    do i=1,nbins
      write(dos_file, *) E(i), dos(i,1), -dos(i,2)
    enddo
  else
    do i=1,nbins
      write(dos_file, *) E(i), dos(i,1)
    enddo
  endif 
  close(dos_file)
end subroutine write_dos


!=============================================================================== 
subroutine setup_energy_scale
!=============================================================================== 
! Sets up all broadening independent DOS concerns
! Calls the relevant dos calculator.
!=============================================================================== 
  implicit none

  integer       :: idos,i,ierr

  allocate(E(1:nbins),stat=ierr)
  if (ierr/=0) call io_error ("cannot allocate E")

  delta_bins=jdos_cutoff_energy/real(nbins-1,dp)
  do idos=1,nbins
     E(idos) = real(idos-1,dp)*delta_bins
  end do

end subroutine setup_energy_scale



!===============================================================================
subroutine allocate_jdos(dos)
!===============================================================================
!===============================================================================
  implicit none

  real(kind=dp), allocatable  :: dos(:,:)
  real(kind=dp), allocatable  :: intdos(:,:)
   
  integer :: ierr 
    
  allocate(dos(nbins,nspins), stat=ierr)
  if(ierr/=0) call io_error("error in allocating dos")
  dos=0.0_dp

 end subroutine allocate_jdos


!===============================================================================
subroutine calculate_jdos(jdos, matrix_weights, weighted_dos)
!===============================================================================

!===============================================================================
 implicit none

 integer :: i,ik,is,ib,idos,ierr,iorb,jb
 integer :: m,n,o,nn
 real(kind=dp) :: dos_temp, cuml, intdos_accum, width
 real(kind=dp) :: grad(1:3), step(1:3), EV(0:4)

 real(kind=dp),intent(out),allocatable, optional    :: weighted_dos(:,:,:)  
 real(kind=dp),intent(in), optional  :: matrix_weights(:,:,:,:)

 real(kind=dp),intent(out),allocatable :: jdos(:,:)
 

 if(linear.or.adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:),dp)/2.0_dp
 if(adaptive) adaptive_smearing=adaptive_smearing*sum(step(:))/3
 if(fixed) width=fixed_smearing
 
 call allocate_jdos(jdos)
 
 do ik=1,nkpoints
    if(iprint>1) then
       if (mod(real(ik,dp),10.0_dp) == 0.0_dp) write(stdout,'(a40,i4,a3,i4,21x,a7)') &
&"Calculating k-point ", ik, " of", nkpoints,"<-- DOS"
    endif
    do is=1,nspins
      do ib=1,vb_max(is)
        do jb=vb_max(is)+1,nbands
           write(stdout,*) ib,jb
           if(linear.or.adaptive) grad(:) = real(band_gradient(jb,jb,:,ik,is)-band_gradient(ib,ib,:,ik,is),dp)*H2ev
           if(linear) call doslin_sub_cell_corners(grad,step,band_energy(jb,is,ik)-band_energy(ib,is,ik),EV)
           if(adaptive) width = sqrt(dot_product(grad,grad))*adaptive_smearing
           
           ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
           ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
           if(finite_bin_correction.and.(width<delta_bins)) width = delta_bins
             
           do idos=1,nbins
             ! The linear method has a special way to calculate the integrated dos
             ! we have to take account for this here.
             if(linear)then
                dos_temp=doslin(EV(0),EV(1),EV(2),EV(3),EV(4),E(idos),cuml)
             else
                dos_temp=gaussian(band_energy(jb,is,ik)-band_energy(ib,is,ik),width,E(idos))&
&*electrons_per_state*kpoint_weight(ik)
              endif

             jdos(idos,is)=jdos(idos,is) + electrons_per_state*kpoint_weight(ik)*dos_temp
                   
          end do
        end do
        end do
     end do
  end do
end subroutine calculate_jdos

end module od_jdos
