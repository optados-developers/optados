! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

!! Progam to calcualte generic spectral
!! R J Nicholls
program gendos
  use od_constants, only: pi, dp
  use od_algorithms, only: gaussian
  use od_io, only: stdout, maxlen, io_file_unit, io_error

  implicit none

  real(kind=dp), allocatable, dimension(:) :: data_energy
  real(kind=dp), allocatable, dimension(:) :: data_value
  real(kind=dp), allocatable, dimension(:, :) :: spectrum   ! Gaussian broadened

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
!  allocate (data(n_transitions, 2), stat=ierr)
  allocate (data_energy(n_transitions), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating data_energy array in gendos')

  allocate (data_value(n_transitions), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating data_value array in gendos')

  data_energy = 0.0_dp
  data_value = 0.0_dp

  do N = 1, n_transitions
    !                          energies    values
    read (gendos_case_unit, *) data_energy(N), data_value(N)
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
  call gaussian_convolute_deltas(data_energy(:), data_value(:), spectrum(:, :), sigma2, gaussian_tol, .false.)

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
  function lorentzian(x0, x, l)
    !! AJ Morris March 2022
    use od_constants, only: dp, pi
    implicit none

    real(kind=dp)             :: lorentzian
    real(dp), intent(in) :: x0, x, l

    lorentzian = l/(pi*(((x0 - x)**2) + (l**2)))
  end function lorentzian

! ==============================================================================
  subroutine lorentzian_convolute_spectrum(spectrum, width)
    !! AJ Morris March 2022

    use od_constants, only: dp
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

! ==============================================================================
  subroutine gaussian_convolute_deltas(input_x, input_y, output, sigma2, gaussian_tol, fast_algor)
    !! AJ Morris March 2022
    use od_algorithms, only: gaussian
    use od_constants, only: dp
    implicit none

    real(dp), dimension(:), intent(in) :: input_x(:) !! x-values (normally positions of delta-fns)
    real(dp), dimension(:), intent(in) :: input_y(:) !! y-values (normally intensities of delta-fns)
    real(dp), dimension(:, :), intent(inout) ::  output(:, :) !! x and y values of smeared spectrum
    real(dp), intent(in) :: sigma2 !! The variance of the Gaussain smearing
    real(dp), intent(in) :: gaussian_tol !! The number of SDs to smear before setting
    !! to zero (for speed)

    logical, optional, intent(in) :: fast_algor !! Use the fast (*2) butterfly algorithm (at the expense of exactitude)

    real(dp) :: centre_of_peak, energy_spacing, minimum_energy, temp_gaussian, nudge
    integer :: centre_bin, start_bin, stop_bin, nbins

    integer :: N, M !! Loop variables

    logical :: fast

    logical :: butterfly !! The butterfly method reduces the number of gaussian calls.
    !! (At the risk of cache thrashing). It is inexact but approaches
    !! the correct answer as nbins --> infity
    !! It assumes that the peak of the gaussian is in the middle of
    !! the bin and there is a small error on the RHS of the gaussian as
    !! a result. You might care -- you might not...

    !! Check for optional argument for the fast method.
    if (present(fast_algor)) then
      butterfly = fast_algor
    else
      butterfly = .false.
    end if

    !! Work out the array sizes for oursleves
    minimum_energy = minval(output(:, 1))

    !! Assume the difference between the first two output energies is the energy spacing.
    energy_spacing = output(2, 1) - output(1, 1)

    nbins = size(output(:, 1))

    !! normal algorthim
    if (.not. butterfly) then
      do N = 1, size(input_x)   !! Loop over each delta function
        centre_of_peak = input_x(N) !! centre of peak
        centre_bin = NINT((input_x(N) - minimum_energy)/energy_spacing)
        start_bin = NINT((input_x(N) - gaussian_tol - minimum_energy)/energy_spacing)
        stop_bin = 2*centre_bin - start_bin

        if (start_bin < 1) start_bin = 1
        if (stop_bin > nbins + 1) stop_bin = nbins + 1

        do M = start_bin, stop_bin   !! Turn each energy value into a function
          output(M, 2) = output(M, 2) + input_y(N)*gaussian(input_x(N), sigma2, output(M, 1))
        end do
      end do

      !! fast algorithm
    elseif (butterfly) then
      do N = 1, size(input_x)   !! Loop over each delta function
        centre_of_peak = input_x(N) !! centre of peak
        centre_bin = NINT((input_x(N) - minimum_energy)/energy_spacing)
        start_bin = NINT((input_x(N) - gaussian_tol - minimum_energy)/energy_spacing)

        !! If we're too near either end of the spectrum then its not worth bothering with the
        !! algorithm.
        fast = .true.
        if (start_bin < 1) Then
          start_bin = 1
          fast = .false.
        end if
        if (stop_bin > nbins + 1) Then
          stop_bin = nbins + 1
          fast = .false.
        end if

        if (fast) then
          do M = start_bin, centre_bin   !! Turn each energy value into a function
            temp_gaussian = gaussian(input_x(N), sigma2, output(M, 1))
            output(M, 2) = output(M, 2) + input_y(N)*temp_gaussian
            output(2*centre_bin - M + 1, 2) = output(2*centre_bin - M + 1, 2) + input_y(N)*temp_gaussian
          end do
        else !! as before
          stop_bin = 2*centre_bin - start_bin
          do M = start_bin, stop_bin   !! Turn each energy value into a function
            output(M, 2) = output(M, 2) + input_y(N)*gaussian(input_x(N), sigma2, output(M, 1))
          end do
        end if
      end do
    end if
  end subroutine gaussian_convolute_deltas

end program gendos
