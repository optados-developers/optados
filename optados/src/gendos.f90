! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
program gendos
  !! Progam to calcualte generic spectral
  !! R J Nicholls
  use od_constants, only: pi, dp
  use od_algorithms, only: gaussian, gaussian_convolute, lorentzian_convolute
  use od_io, only: maxlen, io_file_unit, io_error, io_get_seedname, seedname
  use od_build, only: build_info

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: spectrum_in
  !! The spectrum read from file
  !! (:,1) contains the energies, (:,2) the spectrum
  !! This can be a spectrum to be smeared, or a set of delta functions
  real(kind=dp), allocatable, dimension(:, :) :: spectrum
  !! Intermediate array to hold the spectrum. Same scheme as spectrum_in
  real(kind=dp), allocatable, dimension(:, :) :: spectrum_out
  !! Final array to hold the spectrum. Same scheme as spectrum_in

  real(kind=dp) :: EnergyMinimum
  !! Beginning of spectrum_out scale
  real(kind=dp) :: EnergyMaximum
  !! End of spectrum_out scale
  real(kind=dp) :: EnergySpacing
  !! Energy wdith of each output bin
  real(kind=dp) :: EnergyWindow
  !! When working out EnergyMinimum and EnergyMaximum automoatically this is the
  !! the extra energy added allow and above the highest and lowest delta function
  !! repsectively.
  real(kind=dp) :: GaussianWidth
  !! FWHM of the Gaussian smearing
  real(kind=dp) :: LorentzianWidth
  !! FWHM of the Lorentzian smearing
  real(kind=dp) :: GaussianTolerance
  !! To speed up the calcualtion of the Gaussians, set the Gaussian to zero for
  !! any value GaussianTolerance * HWHM away.
  real(kind=dp) :: LorentzianTolerance
  !! To speed up the calcualtion of the Lorentzian, set the Lorentzian to zero for
  !! any value LorentzianTolerance * HWHM away.

  real(kind=dp) :: sigma
  !! SD of the Gaussian
  real(kind=dp) :: LWidth
  !! Width of Lorentizan (HWHM)

  real(kind=dp) :: gaussian_tol, lorentzian_tol

  integer :: ntransitions, nlines, nbins
  integer :: itransitions, ibin, ilines
  !! Loop variables

  integer :: gendos_input_unit, gendos_output_unit, gendos_case_unit
  integer :: ierr

  character(maxlen) :: cdummy, cdummy2, cdummy3
  character(maxlen) :: data_file

  call io_get_seedname

! If blank set to seedname='--help'
  if (trim(seedname) == '-h' .or. trim(seedname) == '--help') call help_output
  if (trim(seedname) == '-v' .or. trim(seedname) == '--version') call help_output

! Open the input file
  gendos_input_unit = io_file_unit()
  open (unit=gendos_input_unit, action='read', file=trim(seedname)//'.gdi')

  ! Some defaults
  data_file = trim(seedname)//'.dat'
  EnergyMinimum = 1.0_dp    !! Impossible
  EnergyMaximum = -1.0_dp   !! Impossible
  EnergySpacing = 0.01      !! Seems about right (AJM)
  GaussianWidth = 0.0_dp    !! i.e. none
  LorentzianWidth = 0.0_dp  !! i.e none
  GaussianTolerance = 5     !! 5 sigma
  LorentzianTolerance = 10  !! 10x FWHM

  ! Find the number of lines in the input file
  nlines = 0
  do
    read (gendos_input_unit, *, iostat=ierr) cdummy
    if (ierr /= 0) then
      exit
    else
      nlines = nlines + 1
    end if
  end do

  rewind (gendos_input_unit)

  ! Now read the input file
  do ilines = 1, nlines
    read (gendos_input_unit, *, iostat=ierr) cdummy, cdummy2, cdummy3
    select case (trim(cdummy))
    case ("data_file") ! Optional
      read (cdummy3, *) data_file
    case ("ntransitions", "NumTransitions") ! Optional
      read (cdummy3, *) ntransitions
    case ("E_min", "EnergyMin") ! Optional
      read (cdummy3, *) EnergyMinimum
    case ("E_max", "EnergyMax") ! Optional
      read (cdummy3, *) EnergyMaximum
    case ("E_spacing", "EnergySpacing")
      read (cdummy3, *) EnergySpacing
    case ("G_width", "GaussianWidth")
      read (cdummy3, *) GaussianWidth
    case ("L_width", "LorentzianWidth")
      read (cdummy3, *) LorentzianWidth
    case ("G_tol", "GaussianTolerance")
      read (cdummy3, *) GaussianTolerance
    case ("L_tol", "LorentzianTolerance")
      read (cdummy3, *) LorentzianTolerance
    case default
      write (*, *) " Unknown keyword: ", trim(cdummy)
      call io_error('Unknown keyword.')
    end select
  end do

  ! Open the other input data file
  open (unit=gendos_case_unit, action='read', file=trim(data_file))    ! input data file, format energy intensity

  ! Read number of transitions
  ntransitions = 0
  do
    read (gendos_case_unit, *, iostat=ierr) cdummy
    if (ierr /= 0) then
      exit
    else
      ntransitions = ntransitions + 1
    end if
  end do

  rewind (gendos_case_unit)

  ! Read in the data
  allocate (spectrum_in(ntransitions, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum_in array in gendos')

  spectrum_in = 0.0_dp

  do itransitions = 1, ntransitions
    !                          energies    values
    read (gendos_case_unit, *) spectrum_in(itransitions, 1), spectrum_in(itransitions, 2)
  end do

  ! Auto ranging
  if ((EnergyMinimum == 1.0_dp) .and. (EnergyMaximum == -1.0_dp)) then
    ! They haven't been set, let's figure it our for ourselves
    if (GaussianWidth > 0.0) then
      EnergyWindow = 10.0_dp*GaussianWidth
    else ! we did our best
      EnergyWindow = 10.0_dp
    endif
    EnergyMinimum = minval(spectrum_in(:, 1)) - EnergyWindow
    EnergyMaximum = maxval(spectrum_in(:, 1)) + EnergyWindow
  end if

  gendos_output_unit = io_file_unit()
  open (unit=gendos_output_unit, action='write', file=trim(seedname)//'.gdo')

  ! Write out here to catch any defaults that remained unset, or any
  ! vaules we had to work out for ourselves.
  call write_header

  ! Create energy output spectrum
  nbins = NINT((EnergyMaximum - EnergyMinimum)/EnergySpacing)

  allocate (spectrum(nbins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum array in gendos')

  allocate (spectrum_out(nbins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum_out array in gendos')

  spectrum = 0.0_dp
  spectrum_out = 0.0_dp

  ! Put the energies in the array
  do ibin = 1, nbins + 1
    spectrum(ibin, 1) = EnergyMinimum + ((real(ibin) - 1.0_dp)*EnergySpacing)
  end do
  spectrum_out = spectrum

  ! Convert from FWHM notation to sigma^2 notation
  sigma = GaussianWidth/(2*sqrt(2*log(2.0_dp)))

  ! To speed up the gaussian we only do it within tol * sigma of the peak
  gaussian_tol = GaussianTolerance*GaussianWidth*0.5_dp

  ! Smear the spectrum with a guassian of width sigma,
  call gaussian_convolute(spectrum_in(:, :), spectrum(:, :), sigma, gaussian_tol, .false.)

  ! Adds in Lorentzian broadening
  if (LorentzianWidth .gt. 0.0_dp) then
    LWidth = 0.5_dp*LorentzianWidth
    lorentzian_tol = LWidth*LorentzianTolerance
    call lorentzian_convolute(spectrum(:, :), spectrum_out(:, :), LWidth, lorentzian_tol=lorentzian_tol)
  else
    spectrum_out(:, :) = spectrum(:, :)
  end if

  ! Write output
  do ibin = 1, nbins + 1
    write (gendos_output_unit, '(f10.6,f10.6)') spectrum_out(ibin, 1), spectrum_out(ibin, 2)
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
    write (gendos_output_unit, '(a1,a45,a30)') "#", " Input file name :", trim(seedname)//'.gdi'
    write (gendos_output_unit, '(a1,a45,a30)') "#", " Data file name :", trim(data_file)
    write (gendos_output_unit, '(a1,a45,i10)') "#", " Number of transitions in input :", ntransitions
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Minimum energy :", EnergyMinimum
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Maximum energy :", EnergyMaximum
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Energy spacing :", EnergySpacing
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Gaussian FWHM smearing :", GaussianWidth
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Lorentzian FWHM smearing :", LorentzianWidth
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Gaussian precision (in HWHM) :", GaussianTolerance
    write (gendos_output_unit, '(a1,a45,f10.4)') "#", " Lorentzian precision (in HWHM) :", LorentzianTolerance
  end subroutine write_header
end program gendos
