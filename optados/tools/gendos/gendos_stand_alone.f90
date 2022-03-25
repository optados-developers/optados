! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

program gendos

!
! Program to calculate generic spectra
!

  implicit none

! Define the variables
  integer, parameter :: dp = selected_real_kind(15, 300)
  character :: text
!  complex(kind=dp), parameter :: cmplx_i = (0.0_dp,1.0_dp)
  real(kind=dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

  integer :: n_transitions
  integer :: n_bins

  integer :: N
  integer :: M

  real(kind=dp) :: E_min
  real(kind=dp) :: E_max
  real(kind=dp) :: E_spacing
  real(kind=dp) :: G_width
  real(kind=dp) :: l
  real(kind=dp) :: L_width

  real(kind=dp), allocatable, dimension(:, :) :: data
  real(kind=dp), allocatable, dimension(:, :) :: spectrum   ! Gaussian broadened
  real(kind=dp), allocatable, dimension(:) :: spectrum_temp
  real(kind=dp), allocatable, dimension(:, :) :: spectrum_GaL ! Gaussian and Lorentzian broadened
!  real(kind=dp), dimension(3,3) :: recip_cell
!  real(kind=dp), allocatable, dimension(:,:) :: ion_coords

  character(len=30) :: case

! Open the input file
  open (unit=14, action='read', file='gendos.input')
  read (14, *) text, text, case
  read (14, *) text, text, n_transitions
  read (14, *) text, text, E_min
  read (14, *) text, text, E_max
  read (14, *) text, text, E_spacing
  read (14, *) text, text, G_width
  read (14, *) text, text, L_width

! Open the other input and output files
  open (unit=12, action='read', file=case)              ! input data file, format energy intensity
  open (unit=11, action='write', file='gendos.spectra')
!  open(unit=13,action='write',file='')

! Read in the data
  allocate (data(n_transitions, 2))
  data = 0.0_dp
  do N = 1, n_transitions
    read (12, *) data(N, 1), data(N, 2)
  enddo

! Create spectrum (uses gaussian)
  n_bins = NINT((E_max - E_min)/E_spacing)
  allocate (spectrum(n_bins + 1, 2))
  allocate (spectrum_temp(n_bins + 1))
  spectrum = 0.0_dp
  spectrum_temp = 0.0_dp

  Do N = 1, n_bins + 1
    spectrum(N, 1) = E_min + ((N - 1.0_dp)*E_spacing)
  enddo

  do N = 1, n_transitions       ! Loop over energy
    do M = 1, n_bins + 1   ! Turn each energy value into a function
      spectrum_temp(M) = (((4.0_dp*log(2.0_dp))/pi)**(0.5_dp))*(1/G_width)*exp(-4.0_dp*(log(2.0_dp))* &
                                                                               (((spectrum(M, 1) - data(N, 1))/G_width)**2.0_dp))  ! Gaussian
      spectrum(M, 2) = spectrum(M, 2) + (data(N, 2)*spectrum_temp(M))
    end do
  end do                        ! End loop over energy

! Adds in Lorentzian broadening
  if (L_width .gt. 0.00000) then
    allocate (spectrum_GaL(n_bins + 1, 2))
    spectrum_GaL = 0.0_dp
    spectrum_temp = 0.0_dp
    Do N = 1, n_bins + 1
      spectrum_GaL(N, 1) = spectrum(N, 1)
    enddo
    Do N = 1, n_bins + 1       ! Loop over energy
      !         if (E(N_energy) .ge. (LAI_lorentzian_offset + efermi)) then
      !           L_width = 0.5_dp*(LAI_lorentzian_width & ! HWHW of Lorentzian
      !                             + ((E(N_energy) - efermi - LAI_lorentzian_offset)*LAI_lorentzian_scale))
      !         else
      l = 0.5_dp*L_width
      !         end if
      if ((l*pi) .lt. E_spacing) then  ! to get rid of spikes caused by l too small
        l = E_spacing/pi
      endif
      do M = 1, n_bins + 1 ! Turn each energy value into a function
        spectrum_temp(M) = spectrum(N, 2)*l/(pi*(((spectrum(M, 1) - spectrum(N, 1))**2) + (l**2)))*E_spacing  ! Lorentzian
        spectrum_GaL(M, 2) = spectrum_GaL(M, 2) + spectrum_temp(M)
      end do
    end do                        ! End look over energy
  endif

! Write the output file
!  write(11,'(a8,a13)') "# Energy", "Intensity"
  if (L_width .gt. 0.00000) then
    do N = 1, n_bins + 1
      write (11, '(f10.6,f10.6)') spectrum_GaL(N, 1), spectrum_GaL(N, 2)
    enddo
  else
    do N = 1, n_bins + 1
      write (11, '(f10.6,f10.6)') spectrum(N, 1), spectrum(N, 2)
    enddo
  endif

end program gendos
