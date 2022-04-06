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
!=========================================================================!
! Module: algorithms                                                      !
! For low-level algorithms that are not specific to electronic structure  !
!-------------------------------------------------------------------------!
! Modules used:  constants - dp. inv_sqrt_two_pi                          !
!-------------------------------------------------------------------------!
! Key Internal Variables:                                                 !
! Described below                                                         !
!-------------------------------------------------------------------------!
! Necessary conditions:                                                   !
!-------------------------------------------------------------------------!
! Written by Andrew Morris and Chris Picakrd (so far)          11/10/2010 !
!=========================================================================!
module od_algorithms
  use od_constants, only: dp, inv_sqrt_two_pi
  implicit none

  private

  public :: gaussian
  public :: gaussian_convolute
  public :: lorentzian
  public :: lorentzian_convolute
  public :: heap_sort
  public :: utility_lowercase
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: channel_to_am
  public :: algorithms_erf
  public :: algor_dist_array

contains

  function channel_to_am(no)
    implicit none
    character(len=1) :: channel_to_am
    integer, intent(in) :: no

    select case (no)
    case (1)
      channel_to_am = "s"
    case (2)
      channel_to_am = "p"
    case (3)
      channel_to_am = "d"
    case (4)
      channel_to_am = "f"
    case (5)
      channel_to_am = "g"
    end select

  end function channel_to_am

!=========================================================================!
  function gaussian(m, w, x)
!=========================================================================!
! ** Return value of Gaussian(mean=m,width=w) at position x
! I don't know who's this function originally was, CJP? MIJP?
!=========================================================================!
    implicit none

    real(kind=dp), intent(in) :: m, w, x
    real(kind=dp)             :: gaussian

    if (0.5_dp*((x - m)/w)**2 .gt. 30.0_dp) then

      gaussian = 0.0_dp
      return
    else
      gaussian = inv_sqrt_two_pi*exp(-0.5_dp*((x - m)/w)**2)/w
    end if
    return
  end function gaussian

!=========================================================================!
  subroutine heap_sort(num_items, weight)
    !=========================================================================!
    ! This subroutine sorts the list of weights into descending order.        !
    ! On exit, if present, the array of indexes contains the original index   !
    ! of each item.                                                           !
    !                                                                         !
    ! This is a heap sort                                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   num_items (input) :: The number of items to sort                      !
    !   weight (in/out) :: The weights of each item. On exit these are        !
    !                      sorted into descending order.                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Chris Pickard 22nd May 2009                                  !
    !=========================================================================!
    implicit none

    ! Arguments

    integer, intent(in) :: num_items
    real(kind=dp), dimension(num_items), intent(inout) :: weight

    ! Local variables

    integer :: i, ir, j, l ! Loop counters
    real(kind=dp) :: wta

    if (num_items .lt. 2) return

    l = num_items/2 + 1
    ir = num_items

    do
      if (l .gt. 1) then
        l = l - 1
        wta = weight(l)
      else
        wta = weight(ir)
        weight(ir) = weight(1)
        ir = ir - 1
        if (ir .eq. 1) then
          weight(1) = wta
          return
        end if
      end if
      i = l
      j = l + l
20    if (j .le. ir) then
        if (j .lt. ir) then
          if (weight(j) .lt. weight(j + 1)) j = j + 1
        end if
        if (wta .lt. weight(j)) then
          weight(i) = weight(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if
        goto 20
      end if
      weight(i) = wta
    end do

  end subroutine heap_sort

!=========================================================================!
  function utility_lowercase(string)
!=========================================================================!
! Takes a string and converts to lowercase characters
!=========================================================================!

    use od_io, only: maxlen
    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen)        :: utility_lowercase

    integer :: iA, iZ, idiff, ipos, ilett

    iA = ichar('A')
    iZ = ichar('Z')
    idiff = iZ - ichar('z')

    utility_lowercase = string

    do ipos = 1, len(string)
      ilett = ichar(string(ipos:ipos))
      if ((ilett .ge. iA) .and. (ilett .le. iZ)) &
        utility_lowercase(ipos:ipos) = char(ilett - idiff)
    end do

    utility_lowercase = trim(adjustl(utility_lowercase))

    return

  end function utility_lowercase

  !===================================================================
  subroutine utility_frac_to_cart(frac, cart, real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3, 3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i = 1, 3
      cart(i) = real_lat(1, i)*frac(1) + real_lat(2, i)*frac(2) + real_lat(3, i)*frac(3)
    end do

    return

  end subroutine utility_frac_to_cart

  !===================================================================
  subroutine utility_cart_to_frac(cart, frac, recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================
    use od_constants, only: twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3, 3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i = 1, 3
      frac(i) = recip_lat(i, 1)*cart(1) + recip_lat(i, 2)*cart(2) + recip_lat(i, 3)*cart(3)
    end do

    frac = frac/twopi

    return

  end subroutine utility_cart_to_frac

  function algorithms_erf(x)
! Calculate the error function
! From the NSWC Mathematics Subroutine Library
    implicit none
    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: algorithms_erf
    real(kind=dp)             :: C = 0.564189583547756_dp

    real(kind=dp), dimension(1:5) :: A = (/0.771058495001320E-04_dp, -0.133733772997339E-02_dp, &
         & 0.323076579225834E-01_dp, 0.479137145607681E-01_dp, 0.128379167095513E+00_dp/)
    real(kind=dp), dimension(1:3) :: B = (/0.301048631703895E-02_dp, 0.538971687740286E-01_dp, &
         & 0.375795757275549E+00_dp/)
    real(kind=dp), dimension(1:8) :: P = (/-1.36864857382717E-07_dp, 5.64195517478974E-01_dp, &
         & 7.21175825088309E+00_dp, 4.31622272220567E+01_dp, 1.52989285046940E+02_dp, &
         & 3.39320816734344E+02_dp, 4.51918953711873E+02_dp, 3.00459261020162E+02_dp/)
    real(kind=dp), dimension(1:8) :: Q = (/1.00000000000000E+00_dp, 1.27827273196294E+01_dp, &
         & 7.70001529352295E+01_dp, 2.77585444743988E+02_dp, 6.38980264465631E+02_dp, &
         & 9.31354094850610E+02_dp, 7.90950925327898E+02_dp, 3.00459260956983E+02_dp/)
    real(kind=dp), dimension(1:5) :: R = (/2.10144126479064E+00_dp, 2.62370141675169E+01_dp, &
         & 2.13688200555087E+01_dp, 4.65807828718470E+00_dp, 2.82094791773523E-01_dp/)
    real(kind=dp), dimension(1:4) :: S = (/9.41537750555460E+01_dp, 1.87114811799590E+02_dp, &
         & 9.90191814623914E+01_dp, 1.80124575948747E+01_dp/)

    real(kind=dp) :: ax, t, top, bot, x2

    AX = abs(X)
    if (AX <= 0.5_dp) then
      T = X*X
      top = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + 1.0_dp
      bot = ((B(1)*T + B(2))*T + B(3))*T + 1.0_dp
      algorithms_erf = X*(top/bot)
      return
    elseif (AX <= 4.0_dp) then
      top = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX &
           &  + P(6))*AX + P(7))*AX + P(8)
      bot = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX &
           &  + Q(6))*AX + Q(7))*AX + Q(8)
      algorithms_erf = 0.5_dp + (0.5_dp - exp(-X*X)*top/bot)
      if (X < 0.0_dp) algorithms_erf = -algorithms_erf
      return
    elseif (AX <= 5.8_dp) then
      X2 = X*X
      T = 1.0_dp/X2
      top = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
      bot = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + 1.0
      algorithms_erf = (C - top/(X2*bot))/AX
      algorithms_erf = 0.5_dp + (0.5_dp - exp(-X2)*algorithms_erf)
      if (X < 0.0_dp) algorithms_erf = -algorithms_erf
      return
    else
      algorithms_erf = 1.0_dp
      if (X < 0.0_dp) algorithms_erf = -algorithms_erf
    end if

  end function algorithms_erf

  !======================================================
  subroutine algor_dist_array(num_elements, elements_per_node)
    ! Takes the number of elements in an array, num_elements
    ! Returns an array 0,num_nodes-1 which contains the number of
    ! elements that should be on each node.
    ! AJM based on an idea from JRY
    !======================================================
    use od_comms, only: num_nodes
    use od_io, only: io_error
    implicit none
    integer, intent(in)                :: num_elements
    integer, allocatable, intent(out)  :: elements_per_node(:)

    integer :: loop, ierr

    allocate (elements_per_node(0:num_nodes - 1), stat=ierr)
    if (ierr /= 0) call io_error('Error: Problem allocating elements_per_node in algor_dist_array')

    elements_per_node(:) = num_elements/num_nodes

    if (num_elements < num_nodes) then
      call io_error('Fewer kpoints than nodes. Reduce the number of nodes used!')
    end if

    ! Distribute the remainder

    if (elements_per_node(0)*num_nodes .ne. num_elements) then
      do loop = 0, num_elements - elements_per_node(0)*num_nodes - 1
        elements_per_node(loop) = elements_per_node(loop) + 1
      end do
    end if

  end subroutine algor_dist_array

  ! ==============================================================================
  subroutine gaussian_convolute(input, output, sigma2, gaussian_tol, fast_algor)
    !! AJ Morris March 2022
    use od_constants, only: dp
    implicit none

    real(dp), dimension(:, :), intent(in) :: input(:, :) !!
    real(dp), dimension(:, :), intent(out):: output(:, :) !! x and y values of smeared spectrum
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

    do N = 1, size(input, 1)
      write (1, *) input(N, 1), input(N, 2)
    end do

    !! normal algorthim
    if (.not. butterfly) then
      do N = 1, size(input, 1)   !! Loop over each delta function
        centre_of_peak = input(N, 1) !! centre of peak
        centre_bin = NINT((input(N, 1) - minimum_energy)/energy_spacing)
        start_bin = NINT((input(N, 1) - gaussian_tol - minimum_energy)/energy_spacing)
        stop_bin = 2*centre_bin - start_bin

        if (start_bin < 1) start_bin = 1
        if (stop_bin > nbins) stop_bin = nbins
        do M = start_bin, stop_bin   !! Turn each energy value into a function
          output(M, 2) = output(M, 2) + input(N, 2)*gaussian(input(N, 1), sigma2, output(M, 1))
        end do
      end do

      !! fast algorithm
    elseif (butterfly) then
      do N = 1, size(input, 1)   !! Loop over each delta function
        centre_of_peak = input(N, 1) !! centre of peak
        centre_bin = NINT((input(N, 1) - minimum_energy)/energy_spacing)
        start_bin = NINT((input(N, 1) - gaussian_tol - minimum_energy)/energy_spacing)

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
            temp_gaussian = gaussian(input(N, 1), sigma2, output(M, 1))
            output(M, 2) = output(M, 2) + input(N, 2)*temp_gaussian
            output(2*centre_bin - M + 1, 2) = output(2*centre_bin - M + 1, 2) + input(N, 2)*temp_gaussian
          end do
        else !! as before
          stop_bin = 2*centre_bin - start_bin
          do M = start_bin, stop_bin   !! Turn each energy value into a function
            output(M, 2) = output(M, 2) + input(N, 2)*gaussian(input(N, 1), sigma2, output(M, 1))
          end do
        end if
      end do
    end if

    do N = 1, size(output, 1)
      write (3, *) output(N, 1), output(N, 2)
    end do

  end subroutine gaussian_convolute

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
  subroutine lorentzian_convolute(spectrum, width, start, scale)
    !! AJ Morris March 2022

    use od_constants, only: dp, pi
    use od_io, only: io_error
    implicit none

    real(dp) :: energy_spacing

    real(dp), intent(inout), dimension(:, :) :: spectrum
    real(dp), intent(in) :: width !! Halfwidth
    real(dp), allocatable, dimension(:, :) :: spectrum_temp

    real(dp), intent(in), optional :: start, scale

    real(dp) :: l

    integer :: N, M, nbins, ierr

    logical :: energy_dependent_broadening

    if (present(start) .and. present(scale)) then
      energy_dependent_broadening = .true.
    else
      energy_dependent_broadening = .false.
    end if

    nbins = size(spectrum(:, 1))

    ! We need a temporary array to write the output to. We'll deallocate it as
    ! soon as possible.
    allocate (spectrum_temp(1:nbins, 2), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating spectrum_temp array in lorentzian_convolute_spectrum')

    spectrum_temp = 0.0_dp
    spectrum_temp(:, 1) = spectrum(:, 1)

    ! Will want to modify an intent(in) locally. Probably better to do it this way
    ! than end up changing l in the calling function. That would be non-intutive.
    l = width

    energy_spacing = spectrum(2, 1) - spectrum(1, 1)

    if (energy_dependent_broadening) then ! This keeps the if statement out of the loop
      do N = 1, nbins        ! Loop over energy
        !! Annoyingly this equation is for the FULLwidth hence 2.0_dp*l
        if (spectrum(N, 1) .ge. (start)) l = 0.5_dp*(2.0_dp*width + (spectrum(N, 1) - start)*scale)
        if (l*pi .lt. energy_spacing) l = energy_spacing/pi
        do M = 1, nbins ! Turn each energy value into a function
          spectrum_temp(M, 2) = spectrum_temp(M, 2) + &
          &spectrum(N, 2)*energy_spacing*lorentzian(spectrum(M, 1), spectrum(N, 1), l)
        end do
      end do                        ! End look over energy
    else

      !! to get rid of spikes caused by l too small
      if ((l*pi) .lt. energy_spacing) l = energy_spacing/pi

      do N = 1, nbins        ! Loop over energy
        do M = 1, nbins ! Turn each energy value into a function
          !          y  l  E_spacing
          !  -----------------------
          !     pi (x_o - x )^2 + l^2

          spectrum_temp(M, 2) = spectrum_temp(M, 2) + &
          &spectrum(N, 2)*energy_spacing*lorentzian(spectrum(M, 1), spectrum(N, 1), l)
        end do
      end do                        ! End look over energy
    end if

    spectrum = spectrum_temp

    deallocate (spectrum_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating spectrum_GaL array in gendos')

  end subroutine lorentzian_convolute

end module od_algorithms
