! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

!! Progam to calcualte generic spectral
!! R J Nicholls
program gendos
  use od_constants, only: pi, dp
  use od_algorithms, only: gaussian, gaussian_convolute, lorentzian_convolute
  use od_io, only: stdout, maxlen, io_file_unit, io_error, io_get_seedname, seedname
  use od_build, only: build_info

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: input_spectrum, spectrum, spectrum_out   ! Gaussian broadened

  character(maxlen) :: data_file
  real(kind=dp) :: E_min
  real(kind=dp) :: E_max
  real(kind=dp) :: E_spacing
  real(kind=dp) :: G_width
  real(kind=dp) :: L_width
  real(kind=dp) :: tol

  real(kind=dp) :: l

  real(kind=dp) :: sigma2
  real(kind=dp) :: centre_of_peak
  real(kind=dp) :: gaussian_tol
  real(kind=dp) :: energy_spacing, minimum_energy

  integer :: n_transitions, n_lines, ilines
  integer :: n_bins
  integer :: N
  integer :: M
  integer :: gendos_input_unit
  integer :: gendos_output_unit
  integer :: gendos_case_unit
  integer :: ierr, iostat
  integer :: Start_bin, stop_bin, centre_bin

  character(maxlen) :: cdummy, cdummy2, cdummy3

  call io_get_seedname

! If blank set to seedname='--help'
  if (trim(seedname) == '-h' .or. trim(seedname) == '--help') call help_output
  if (trim(seedname) == '-v' .or. trim(seedname) == '--version') call help_output

! Open the input file
  gendos_input_unit = io_file_unit()
  open (unit=gendos_input_unit, action='read', file=trim(seedname)//'.gdi')

  ! Some defaults
  data_file = trim(seedname)//'.dat'
  E_min = 1.0_dp  !! Impossible
  E_max = -1.0_dp !! Impossible
  E_spacing = 0.01    !! Seems about right (AJM)
  G_width = 0.0_dp  !! i.e. none
  L_width = 0.0_dp  !! i.e none
  tol = 5       !! 5 sigma

  n_lines = 0
  do
    read (gendos_input_unit, *, iostat=ierr) cdummy
    if (ierr /= 0) then
      exit
    else
      n_lines = n_lines + 1
    end if
  end do

  rewind (gendos_input_unit)

  do ilines = 1, n_lines
    read (gendos_input_unit, *, iostat=ierr) cdummy, cdummy2, cdummy3
    select case (trim(cdummy))
    case ("data_file") ! Optional
      read (cdummy3, *) data_file
    case ("n_transitions", "NumTransitions") ! Optional
      read (cdummy3, *) n_transitions
    case ("E_min", "EnergyMin") ! Optional
      read (cdummy3, *) E_min
    case ("E_max", "EnergyMax") ! Optional
      read (cdummy3, *) E_max
    case ("E_spacing", "EnergySpacing")
      read (cdummy3, *) E_spacing
    case ("G_width", "GaussianWidth")
      read (cdummy3, *) G_width
    case ("L_width", "LorentzianWidth")
      read (cdummy3, *) L_width
    case ("Tol", "GaussianTolerance")
      read (cdummy3, *) tol
    case default
      write (*, *) " Unknown keyword: ", trim(cdummy)
      call io_error('Unknown keyword.')
    end select
  end do

  !! Open the other input and output files
  open (unit=gendos_case_unit, action='read', file=trim(data_file))    ! input data file, format energy intensity

  ! Read number of transitions
  n_transitions = 0
  do
    read (gendos_case_unit, *, iostat=ierr) cdummy
    if (ierr /= 0) then
      exit
    else
      n_transitions = n_transitions + 1
    end if
  end do

  rewind (gendos_case_unit)

  ! Read in the data
  allocate (input_spectrum(n_transitions, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating input_spectrum array in gendos')

  input_spectrum = 0.0_dp

  do N = 1, n_transitions
    !                          energies    values
    read (gendos_case_unit, *) input_spectrum(N, 1), input_spectrum(N, 2)
  end do

  ! Auto ranging
  if ((E_min == 1.0_dp) .and. (E_max == -1.0_dp)) then
    ! They haven't been set, let's figure it our for ourselves
    E_min = minval(input_spectrum(:, 1)) - 10.0_dp*G_width
    E_max = maxval(input_spectrum(:, 1)) + 10.0_dp*G_width
  end if

  gendos_output_unit = io_file_unit()
  open (unit=gendos_output_unit, action='write', file=trim(seedname)//'.gdo')

  !! Write out here to catch defaults
  call write_header

  ! Create spectrum
  n_bins = NINT((E_max - E_min)/E_spacing)

  allocate (spectrum(n_bins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum array in gendos')
  allocate (spectrum_out(n_bins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum_out array in gendos')

  spectrum = 0.0_dp
  spectrum_out = 0.0_dp

  ! Put the energies in the array
  do N = 1, n_bins + 1
    spectrum(N, 1) = E_min + ((N - 1.0_dp)*E_spacing)
  end do
  spectrum_out = spectrum

  ! Convert from FWHM notation to sigma^2 notation
  sigma2 = G_width/(2*sqrt(2*log(2.0_dp)))

  ! To speed up the gaussian we only do it within tol * sigma of the peak
  gaussian_tol = tol*sqrt(sigma2)

  ! Smear the spectrum with a guassian of width sigma2,
  call gaussian_convolute(input_spectrum(:, :), spectrum(:, :), sigma2, gaussian_tol, .false.)

  ! Adds in Lorentzian broadening
  if (L_width .gt. 0.0_dp) then
    l = 0.5_dp*L_width
    call lorentzian_convolute(spectrum(:, :), spectrum_out(:, :), l)
  end if

  do N = 1, n_bins + 1
    write (gendos_output_unit, '(f10.6,f10.6)') spectrum_out(N, 1), spectrum_out(N, 2)
  end do

contains
  subroutine help_output
    use od_constants, only: optados_version, copyright
    implicit none
    write (*, *)
    write (*, *) " Gendos version ", trim(build_info%build)
    write (*, *)
    write (*, *) " Andrew J. Morris, R. J. Nicholls, C. J. Pickard and J. R. Yates", trim(copyright)
    write (*, *) " Usage: gendos <seedname>"
    write (*, *)
    stop
  end subroutine help_output

  subroutine write_header
    implicit none
    write (gendos_output_unit, '(a78)') "# +==========================================================================+"
    write (gendos_output_unit, '(a78)') "# |                                                                          | "
    write (gendos_output_unit, '(a78)') "# |             GGG    EEEE  N   N   DDD    OOO    SSS                       | "
    write (gendos_output_unit, '(a78)') "# |            G       E     NN  N   D  D  O   O  S                          | "
    write (gendos_output_unit, '(a78)') "# |            G  GGG  EEE   N N N   D  D  O   O   SS                        | "
    write (gendos_output_unit, '(a78)') "# |            G   G   E     N  NN   D  D  O   O     S                       | "
    write (gendos_output_unit, '(a78)') "# |             GGG    EEEE  N   N   DDD    OOO   SSS                        | "
    write (gendos_output_unit, '(a78)') "# |                                                                          | "
    write (gendos_output_unit, '(a78)') "# +--------------------------------------------------------------------------+ "
    write (gendos_output_unit, '(a1,a40,a30)') "#", " Input file name :", trim(seedname)//'.gdi'
    write (gendos_output_unit, '(a1,a40,a30)') "#", " Data file name :", trim(data_file)
    write (gendos_output_unit, '(a1,a40,i10)') "#", " Number of transitions in input :", n_transitions
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Minimum energy :", E_min
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Maximum energy :", E_max
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Energy spacing :", E_spacing
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Gaussian FWHM smearing :", G_width
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Lorentzian FWHM smearing :", L_width
    write (gendos_output_unit, '(a1,a40,f10.4)') "#", " Gaussian precision (in SDs) :", Tol
  end subroutine write_header
end program gendos
