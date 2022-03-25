! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

!! Progam to calcualte generic spectral
!! R J Nicholls
program gendos
  use od_constants, only: pi, dp
  use od_algorithms, only: gaussian
  use od_io, only: stdout, maxlen, io_file_unit, io_error

  implicit none

  real(kind=dp), allocatable, dimension(:, :) :: data
  real(kind=dp), allocatable, dimension(:, :) :: spectrum   ! Gaussian broadened
  real(kind=dp), allocatable, dimension(:, :) :: spectrum_GaL ! Gaussian and Lorentzian broadened

  real(kind=dp) :: E_min
  real(kind=dp) :: E_max
  real(kind=dp) :: E_spacing
  real(kind=dp) :: G_width
  real(kind=dp) :: l
  real(kind=dp) :: L_width
  real(kind=dp) :: spectrum_temp
  real(kind=dp) :: sigma_2

  integer :: n_transitions
  integer :: n_bins
  integer :: N
  integer :: M
  integer :: gendos_input_unit
  integer :: gendos_output_unit
  integer :: gendos_case_unit
  integer :: ierr

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
  allocate (data(n_transitions, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating data array in gendos')

  data = 0.0_dp

  do N = 1, n_transitions
    read (gendos_case_unit, *) data(N, 1), data(N, 2)
  end do

  ! Create spectrum (uses gaussian)
  n_bins = NINT((E_max - E_min)/E_spacing)

  allocate (spectrum(n_bins + 1, 2), stat=ierr)
  if (ierr /= 0) call io_error('Error allocating spectrum array in gendos')

  spectrum = 0.0_dp
  spectrum_temp = 0.0_dp

  ! Put the energies in the  array
  do N = 1, n_bins + 1
    spectrum(N, 1) = E_min + ((N - 1.0_dp)*E_spacing)
  end do

  ! Convert from FWHM notation to sigma^2 notation
  sigma2 = G_width/(2*sqrt(2*log(2.0_dp)))

  do N = 1, n_transitions       ! Loop over energy
    do M = 1, n_bins + 1   ! Turn each energy value into a function

      !function gaussian(m, w, x)
      !gaussian = inv_sqrt_two_pi*exp(-0.5_dp*((x - m)/w)**2)/w

      ! need to convert FWHM to variance
      !  m = data(N, 1)
      !  w = sigma_2
      !  x = spectrum(M, 1)

      !       spectrum_temp = gaussian(data(N, 1), sigma2, spectrum(M, 1))
      !     __________
      !    / 4*log(2)       1         {            (    S - D   )**2  }
      !   / ---------  * -------  exp { -4 log(2)  (------------)     }
      ! \/      Pi       G_width      {            ( G_width    )     }

      !   spectrum_temp = (((4.0_dp*log(2.0_dp))/pi)**(0.5_dp))*(1/G_width)*exp(-4.0_dp*(log(2.0_dp))* &
      !            & (((spectrum(M, 1) - data(N, 1))/G_width)**2.0_dp))  ! Gaussian

      ! Option to truncate gaussian
      spectrum(M, 2) = spectrum(M, 2) + data(N, 2)*gaussian(data(N, 1), sigma2, spectrum(M, 1))
    end do
  end do

! Adds in Lorentzian broadening
  if (L_width .gt. 0.0_dp) then
    allocate (spectrum_GaL(n_bins + 1, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating spectrum_GaL array in gendos')

    spectrum_GaL = 0.0_dp
    spectrum_temp = 0.0_dp

    spectrum_GaL(:, 1) = spectrum(:, 1)

    l = 0.5_dp*L_width

    do N = 1, n_bins + 1       ! Loop over energy
      if ((l*pi) .lt. E_spacing) l = E_spacing/pi  ! to get rid of spikes caused by l too small

      do M = 1, n_bins + 1 ! Turn each energy value into a function
        spectrum_temp = spectrum(N, 2)*l/(pi*(((spectrum(M, 1) - spectrum(N, 1))**2) + (l**2)))*E_spacing  ! Lorentzian
        spectrum_GaL(M, 2) = spectrum_GaL(M, 2) + spectrum_temp
      end do
    end do                        ! End look over energy
  end if

  if (L_width .gt. 0.00_dp) then
    do N = 1, n_bins + 1
      write (gendos_output_unit, '(f10.6,f10.6)') spectrum_GaL(N, 1), spectrum_GaL(N, 2)
    end do
  else
    do N = 1, n_bins + 1
      write (gendos_output_unit, '(f10.6,f10.6)') spectrum(N, 1), spectrum(N, 2)
    end do
  end if

end program gendos
