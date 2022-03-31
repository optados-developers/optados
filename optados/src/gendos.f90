! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

!! Progam to calcualte generic spectral
!! R J Nicholls
program gendos
  use od_constants, only: pi, dp
  use od_algorithms, only: gaussian, gaussian_convolute
  use od_io, only: stdout, maxlen, io_file_unit, io_error

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: input_spectrum, spectrum, spectrum_tmp   ! Gaussian broadened

  real(kind=dp) :: E_min
  real(kind=dp) :: E_max
  real(kind=dp) :: E_spacing
  real(kind=dp) :: G_width
  real(kind=dp) :: l
  real(kind=dp) :: L_width
  real(kind=dp) :: sigma2
  real(kind=dp) :: centre_of_peak
  real(kind=dp) :: gaussian_tol
  real(kind=dp) :: energy_spacing, minimum_energy

  integer :: n_transitions
  integer :: n_bins
  integer :: N
  integer :: M
  integer :: gendos_input_unit
  integer :: gendos_output_unit
  integer :: gendos_case_unit
  integer :: ierr
  integer :: Start_bin, stop_bin, centre_bin

  character(maxlen) :: cdummy
  character(maxlen) :: case

  gendos_input_unit = io_file_unit()

! Open the input file
  open (unit=gendos_input_unit, action='read', file='gendos.input')
  read (gendos_input_unit, *) cdummy, cdummy, case
  read (gendos_input_unit, *) cdummy, cdummy, n_transitions
  read (gendos_input_unit, *) cdummy, cdummy, E_min
  read (gendos_input_unit, *) cdummy, cdummy, E_max
  read (gendos_input_unit, *) cdummy, cdummy, E_spacing
  read (gendos_input_unit, *) cdummy, cdummy, G_width
  read (gendos_input_unit, *) cdummy, cdummy, L_width

  gendos_output_unit = io_file_unit()

! Open the other input and output files
  open (unit=gendos_case_unit, action='read', file=case)      ! input data file, format energy intensity
  open (unit=gendos_output_unit, action='write', file='gendos.spectra')

! Read in the data
  allocate (input_spectrum(n_transitions, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating input_spectrum array in gendos')

  input_spectrum = 0.0_dp

  do N = 1, n_transitions
    !                          energies    values
    read (gendos_case_unit, *) input_spectrum(N, 1), input_spectrum(N, 2)
  end do

  ! Create spectrum
  n_bins = NINT((E_max - E_min)/E_spacing)

  allocate (spectrum(n_bins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum array in gendos')

  spectrum = 0.0_dp

  ! Put the energies in the  array
  do N = 1, n_bins + 1
    spectrum(N, 1) = E_min + ((N - 1.0_dp)*E_spacing)
  end do

  ! Convert from FWHM notation to sigma^2 notation
  sigma2 = G_width/(2*sqrt(2*log(2.0_dp)))

  ! To speed up the gaussian we only do it within 2 * sigma of the peak
  gaussian_tol = 2*sqrt(sigma2)

  ! Smear the spectrum with a guassian of width sigma2,
  call gaussian_convolute(input_spectrum(:, :), spectrum(:, :), sigma2, gaussian_tol, .false.)

! Adds in Lorentzian broadening
  if (L_width .gt. 0.0_dp) then
    l = 0.5_dp*L_width
    call lorentzian_convolute_spectrum(spectrum(:, :), l)
  end if

  do N = 1, n_bins + 1
    write (gendos_output_unit, '(f10.6,f10.6)') spectrum(N, 1), spectrum(N, 2)
  end do

contains

! ==============================================================================
  subroutine lorentzian_convolute_spectrum(spectrum, width)
    !! AJ Morris March 2022

    use od_constants, only: dp
    use od_algorithms, only: lorentzian
    implicit none

    real(dp) :: energy_spacing

    real(dp), intent(inout), dimension(:, :) :: spectrum
    real(dp), intent(in) :: width
    real(dp), allocatable, dimension(:, :) :: spectrum_temp

    real(dp) :: l

    integer :: N, M, nbins

    nbins = size(spectrum(:, 1))

    ! We need a temporary array to write the output to. We'll deallocate it as
    ! soon as possible.
    allocate (spectrum_temp(n_bins + 1, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating spectrum_temp array in lorentzian_convolute_spectrum')

    spectrum_temp = 0.0_dp
    spectrum_temp(:, 1) = spectrum(:, 1)

    ! Will want to modify an intent(in) locally. Probably better to do it this way
    ! than end up changing l in the calling function. That would be non-intutive.
    l = width

    energy_spacing = spectrum(2, 1) - spectrum(1, 1)

    !! to get rid of spikes caused by l too small
    if ((l*pi) .lt. energy_spacing) l = energy_spacing/pi

    do N = 1, nbins + 1       ! Loop over energy
      do M = 1, nbins + 1 ! Turn each energy value into a function
        !
        !          y  l  E_spacing
        !  -----------------------
        !     pi (x_o - x )^2 + l^2

        spectrum_temp(M, 2) = spectrum_temp(M, 2) + spectrum(N, 2)*energy_spacing*lorentzian(spectrum(M, 1), spectrum(N, 1), l)
      end do
    end do                        ! End look over energy

    spectrum = spectrum_temp

    deallocate (spectrum_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating spectrum_GaL array in gendos')

  end subroutine lorentzian_convolute_spectrum

end program gendos
