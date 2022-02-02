!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
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
module xmgrace_utils
  ! use od_constants, only : dp

  implicit none

  public :: xmgu_setup
  public :: xmgu_legend
  public :: xmgu_title
  public :: xmgu_subtitle
  public :: xmgu_axis
  public :: xmgu_data
  public :: xmgu_data_header
  public :: xmgu_vertical_line  !(unit,x_coord,y_max,y_min)

  private

contains

  !==================================================================!
  subroutine xmgu_setup(unit)
    !==================================================================!
    use od_io, only: io_date
    implicit none
    integer, intent(in)      :: unit
    character(len=9) :: ctime             ! Temp. time string
    character(len=11):: cdate             ! Temp. date string

    call io_date(cdate, ctime)
    write (unit, *) '# Grace project file'
    write (unit, *) '# Autogenrated by OptaDOS on ', cdate, ' at ', ctime
    write (unit, *) '@version 50122'
    write (unit, *) '@page size 792, 612'
    write (unit, *) '@page scroll 5%'
    write (unit, *) '@page inout 5%'
    write (unit, *) '@link page off'

    write (unit, *) '@map font 0 to "Times-Roman", "Times-Roman"'
    write (unit, *) '@map font 1 to "Times-Italic", "Times-Italic"'
    write (unit, *) '@map font 2 to "Times-Bold", "Times-Bold"'
    write (unit, *) '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"'
    write (unit, *) '@map font 4 to "Helvetica", "Helvetica"'
    write (unit, *) '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"'
    write (unit, *) '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"'
    write (unit, *) '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"'
    write (unit, *) '@map font 8 to "Courier", "Courier"'
    write (unit, *) '@map font 9 to "Courier-Oblique", "Courier-Oblique"'
    write (unit, *) '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write (unit, *) '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"'
    write (unit, *) '@map font 12 to "Symbol", "Symbol"'
    write (unit, *) '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
    write (unit, *) '@map color 0 to (255, 255, 255), "white"'
    write (unit, *) '@map color 1 to (0, 0, 0), "black"'
    write (unit, *) '@map color 2 to (255, 0, 0), "red"'
    write (unit, *) '@map color 3 to (0, 255, 0), "green"'
    write (unit, *) '@map color 4 to (0, 0, 255), "blue"'
    write (unit, *) '@map color 5 to (255, 255, 0), "yellow"'
    write (unit, *) '@map color 6 to (188, 143, 143), "brown"'
    write (unit, *) '@map color 7 to (220, 220, 220), "grey"'
    write (unit, *) '@map color 8 to (148, 0, 211), "violet"'
    write (unit, *) '@map color 9 to (0, 255, 255), "cyan"'
    write (unit, *) '@map color 10 to (255, 0, 255), "magenta"'
    write (unit, *) '@map color 11 to (255, 165, 0), "orange"'
    write (unit, *) '@map color 12 to (114, 33, 188), "indigo"'
    write (unit, *) '@map color 13 to (103, 7, 72), "maroon"'
    write (unit, *) '@map color 14 to (64, 224, 208), "turquoise"'

    write (unit, *) '@default linewidth 1.0'
    write (unit, *) '@default linestyle 1'
    write (unit, *) '@default color 1'
    write (unit, *) '@default pattern 1'
    write (unit, *) '@default font 0'
    write (unit, *) '@default char size 1.000000'
    write (unit, *) '@default symbol size 1.000000'
    write (unit, *) '@default sformat "%.8g"'

    write (unit, *) '@background color 0'
    write (unit, *) '@page background fill on'

    write (unit, *) '@g0 on'
    write (unit, *) '@g0 hidden false'
    write (unit, *) '@g0 type XY'
    write (unit, *) '@g0 stacked false'
    write (unit, *) '@g0 bar hgap 0.000000'
    write (unit, *) '@g0 fixedpoint off'
    write (unit, *) '@g0 fixedpoint type 0'
    write (unit, *) '@g0 fixedpoint xy 0.000000, 0.000000'
    write (unit, *) '@g0 fixedpoint format general general'
    write (unit, *) '@g0 fixedpoint prec 6, 6'

!     write(unit,*) 'autoscale'
!     write(unit,*) 'world xmin 0.0'
!     write(unit,*) 'world xmax 1.0'
    !    write(unit,*) 'world ymax 0.1'

    !    write(unit,*) 'frame linestyle 1'
    !   write(unit,*) 'frame linewidth 2.5'
    !   write(unit,*) 'frame color 1'
    !   write(unit,*) 'frame pattern 1'
    !   write(unit,*) 'frame background color 0'
    !   write(unit,*) 'frame background pattern 0'

    !   write(unit,*) 'altxaxis  off'
    !   write(unit,*) 'altyaxis  off'

  end subroutine xmgu_setup
  !==================================================================!

  !==================================================================!
  subroutine xmgu_legend(unit)
    !==================================================================!
    implicit none
    integer, intent(in)      :: unit

    write (unit, *) '@    legend on'
    write (unit, *) '@    legend loctype view'
    write (unit, *) '@    legend 0.85, 0.8'
    write (unit, *) '@    legend box color 1'
    write (unit, *) '@    legend box pattern 1'
    write (unit, *) '@    legend box linewidth 2.0'
    write (unit, *) '@    legend box linestyle 1'
    write (unit, *) '@    legend box fill color 0'
    write (unit, *) '@    legend box fill pattern 1'
    write (unit, *) '@    legend font 4'
    write (unit, *) '@    legend char size 1.000000'
    write (unit, *) '@    legend color 1'
    write (unit, *) '@    legend length 4'
    write (unit, *) '@    legend vgap 1'
    write (unit, *) '@    legend hgap 1'
    write (unit, *) '@    legend invert false'

  end subroutine xmgu_legend
  !==================================================================!

  !==================================================================!
  subroutine xmgu_title(unit, min_x, max_x, min_y, max_y, title)
    !==================================================================!
    use od_constants, only: dp
    implicit none

    integer, intent(in)      :: unit
    character(*), intent(in) :: title

    real(dp), intent(in) :: min_x, max_x, min_y, max_y

    character(20) :: min_x_char, min_y_char, max_x_char, max_y_char

    write (min_x_char, '(F20.6)') min_x
    write (max_x_char, '(F20.6)') max_x
    write (min_y_char, '(F20.6)') min_y
    write (max_y_char, '(F20.6)') max_y

    min_x_char = trim(adjustl((min_x_char)))
    max_x_char = trim(adjustl((max_x_char)))
    min_y_char = trim(adjustl((min_y_char)))
    max_y_char = trim(adjustl((max_y_char)))

    write (unit, *) '    @with g0'
    write (unit, *) '@    world '//trim(min_x_char)//', '//trim(min_y_char)&
     &//', '//trim(max_x_char)//', '//trim(max_y_char)
    write (unit, *) '@    stack world 0, 0, 0, 0'
    write (unit, *) '@    znorm 1'
    write (unit, *) '@    view 0.150000, 0.150000, 1.150000, 0.850000'
    write (unit, *) '@    title "'//trim(title)//'"'
    write (unit, *) '@    title font 4'
    write (unit, *) '@    title size 1.500000'
    write (unit, *) '@    title color 1'

    ! write(unit,*) 'title "'//trim(title)//'"'
    ! write(unit,*) "title font 4"
    ! write(unit,*) "title size 1.500000"
    ! write(unit,*) "title color 1"

  end subroutine xmgu_title
  !==================================================================!

  !==================================================================!
  subroutine xmgu_subtitle(unit, subtitle)
    !==================================================================!
    implicit none

    integer, intent(in)      :: unit
    character(*), intent(in) :: subtitle

    write (unit, *) '@    subtitle "'//trim(subtitle)//'"'
    write (unit, *) '@    subtitle font 4'
    write (unit, *) '@    subtitle size 1.000000'
    write (unit, *) '@    subtitle color 1'

    !   write(unit,*) 'subtitle '//trim(subtitle)//'"'
    !   write(unit,*) "subtitle font 4"
    !   write(unit,*) "subtitle size 1.000000"
    !   write(unit,*) "subtitle color 1"

  end subroutine xmgu_subtitle
  !==================================================================!

  !==================================================================!
  subroutine xmgu_axis(unit, axis, label)
    !==================================================================!
    implicit none

    integer, intent(in)      :: unit
    character(*), intent(in) :: label
    character(*), intent(in)  :: axis

    character(5) :: axis_name

    if (axis == "x") then
      axis_name = "xaxis"
    elseif (axis == "y") then
      axis_name = "yaxis"
    end if

    write (unit, *) '@    '//trim(axis_name)//' on'
    write (unit, *) '@    '//trim(axis_name)//' type zero false'
    write (unit, *) '@    '//trim(axis_name)//' offset 0.000000 , 0.000000'
    write (unit, *) '@    '//trim(axis_name)//' bar on'
    write (unit, *) '@    '//trim(axis_name)//' bar color 1'
    write (unit, *) '@    '//trim(axis_name)//' bar linestyle 1'
    write (unit, *) '@    '//trim(axis_name)//' bar linewidth 2.5'
    write (unit, *) '@    '//trim(axis_name)//' label "'//trim(label)//'"'
    write (unit, *) '@    '//trim(axis_name)//' label layout para'
    write (unit, *) '@    '//trim(axis_name)//' label place auto'
    write (unit, *) '@    '//trim(axis_name)//' label char size 1.000000'
    write (unit, *) '@    '//trim(axis_name)//' label font 4'
    write (unit, *) '@    '//trim(axis_name)//' label color 1'
    write (unit, *) '@    '//trim(axis_name)//' label place normal'
    write (unit, *) '@    '//trim(axis_name)//' tick on'
    write (unit, *) '@ autoticks'
    !  write(unit,*)'@    '//trim(axis_name)//' tick major 10'
    !  write(unit,*)'@    '//trim(axis_name)//' tick minor ticks 1'
    !  write(unit,*)'@    '//trim(axis_name)//' tick default 6'
    write (unit, *) '@    '//trim(axis_name)//' tick place rounded true'
    write (unit, *) '@    '//trim(axis_name)//' tick in'
    write (unit, *) '@    '//trim(axis_name)//' tick major size 1.000000'
    write (unit, *) '@    '//trim(axis_name)//' tick major color 1'
    write (unit, *) '@    '//trim(axis_name)//' tick major linewidth 2.5'
    write (unit, *) '@    '//trim(axis_name)//' tick major linestyle 1'
    write (unit, *) '@    '//trim(axis_name)//' tick major grid off'
    write (unit, *) '@    '//trim(axis_name)//' tick minor color 1'
    write (unit, *) '@    '//trim(axis_name)//' tick minor linewidth 2.5'
    write (unit, *) '@    '//trim(axis_name)//' tick minor linestyle 1'
    write (unit, *) '@    '//trim(axis_name)//' tick minor grid off'
    write (unit, *) '@    '//trim(axis_name)//' tick minor size 0.500000'
    write (unit, *) '@    '//trim(axis_name)//' ticklabel on'
    write (unit, *) '@    '//trim(axis_name)//' ticklabel format general'
    write (unit, *) '@    '//trim(axis_name)//' ticklabel font 4'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel prec 5'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel formula ""'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel append ""'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel prepend ""'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel angle 0'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel skip 0'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel stagger 0'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel place normal'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel offset auto'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel offset 0.000000 , 0.010000'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel start type auto'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel start 0.000000'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel stop type auto'
    ! 1write(unit,*)'@    '//trim(axis_name)//' ticklabel stop 0.000000'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel char size 1.000000'
    ! write(unit,*)'@    '//trim(axis_name)//' ticklabel color 1'
    ! write(unit,*)'@    '//trim(axis_name)//' tick place both'
    ! write(unit,*)'@    '//trim(axis_name)//' tick spec type none'
!    write(unit,*) trim(axis_name)//" type zero true"
!    write(unit,*) trim(axis_name)//" bar on"
!    write(unit,*) trim(axis_name)//" bar color 1"
!    write(unit,*) trim(axis_name)//" bar linestyle 1"
!    write(unit,*) trim(axis_name)//" bar linewidth 2.5"
!    write(unit,*) trim(axis_name)//' label "'//trim(label)//'"'
!    write(unit,*) trim(axis_name)//" label layout para"
!    write(unit,*) trim(axis_name)//" label place auto"
!    write(unit,*) trim(axis_name)//" label char size 1.000000"
!    write(unit,*) trim(axis_name)//" label font 4"
!    write(unit,*) trim(axis_name)//" label color 1"
!    write(unit,*) trim(axis_name)//" label place normal"
!    write(unit,*) trim(axis_name)//" tick on"
!    write(unit,*) trim(axis_name)//" tick major 0.2"
!    write(unit,*) trim(axis_name)//" tick minor ticks 1"
!    write(unit,*) trim(axis_name)//" tick default 6"
!    write(unit,*) trim(axis_name)//" tick place rounded true"
!    write(unit,*) trim(axis_name)//" tick in"
!    write(unit,*) trim(axis_name)//" tick major size 1.000000"
!    write(unit,*) trim(axis_name)//" tick major color 1"
!    write(unit,*) trim(axis_name)//" tick major linewidth 2.5"
!    write(unit,*) trim(axis_name)//" tick major linestyle 1"
!    write(unit,*) trim(axis_name)//" tick major grid off"
!    write(unit,*) trim(axis_name)//" tick minor color 1"
!    write(unit,*) trim(axis_name)//" tick minor linewidth 2.5"
!    write(unit,*) trim(axis_name)//" tick minor linestyle 1"
!    write(unit,*) trim(axis_name)//" tick minor grid off"
!    write(unit,*) trim(axis_name)//" tick minor size 0.500000"
!    write(unit,*) trim(axis_name)//" ticklabel on"
!    write(unit,*) trim(axis_name)//" ticklabel format general"
!    write(unit,*) trim(axis_name)//" ticklabel prec 5"
!    write(unit,*) trim(axis_name)//" ticklabel start type spec"
!    write(unit,*) trim(axis_name)//" ticklabel start 0.200000"
!    write(unit,*) trim(axis_name)//" ticklabel stop type spec"
    !   write(unit,*) trim(axis_name)//" ticklabel stop 0.800000"
!    write(unit,*) trim(axis_name)//" ticklabel font 4"

  end subroutine xmgu_axis
  !==================================================================!

  !==================================================================!
  subroutine xmgu_data(unit, field, x_data, y_data)
    !==================================================================!

    use od_constants, only: dp
    use od_io, only: io_error
    implicit none
    integer, intent(in)      :: unit
    integer, intent(in)      :: field

    character(4)              :: char_field
    real(dp), intent(in)       :: x_data(:)
    real(dp), intent(in)       :: y_data(:)
    integer :: i

    write (char_field, '(I4)') field

    char_field = trim("s"//adjustl(char_field))

    write (unit, *) '@target G0.'//trim(char_field)
    write (unit, *) '@type xy'

    if (size(x_data, 1) .ne. size(y_data, 1)) call io_error("xmgu_data: x and y axes have different number of elements")

    do i = 1, size(x_data, 1)
      write (unit, *) x_data(i), y_data(i)
    end do
    write (unit, *) "&"

  end subroutine xmgu_data

  !==================================================================!
  subroutine xmgu_data_header(unit, field, colour, legend)
    !==================================================================!
    integer, intent(in)      :: unit
    integer, intent(in)      :: field, colour
    character(4)              :: char_field, char_colour
    character(*), intent(in)  :: legend

    write (char_field, '(I4)') field
    write (char_colour, '(I4)') colour

    char_field = trim("s"//adjustl(char_field))

    write (unit, *) '@    '//trim(char_field)//' hidden false'
    write (unit, *) '@    '//trim(char_field)//' type xy'
    write (unit, *) '@    '//trim(char_field)//' symbol 0'
    write (unit, *) '@    '//trim(char_field)//' symbol size 1.000000'
    write (unit, *) '@    '//trim(char_field)//' symbol color '//trim(char_colour)
    write (unit, *) '@    '//trim(char_field)//' symbol pattern 1'
    write (unit, *) '@    '//trim(char_field)//' symbol fill color 2'
    write (unit, *) '@    '//trim(char_field)//' symbol fill pattern 0'
    write (unit, *) '@    '//trim(char_field)//' symbol linewidth 1.0'
    write (unit, *) '@    '//trim(char_field)//' symbol linestyle 1'
    write (unit, *) '@    '//trim(char_field)//' symbol char 65'
    write (unit, *) '@    '//trim(char_field)//' symbol char font 0'
    write (unit, *) '@    '//trim(char_field)//' symbol skip 0'
    write (unit, *) '@    '//trim(char_field)//' line type 1'
    write (unit, *) '@    '//trim(char_field)//' line linestyle 1'
    write (unit, *) '@    '//trim(char_field)//' line linewidth 2.0'
    write (unit, *) '@    '//trim(char_field)//' line color  '//trim(char_colour)
    write (unit, *) '@    '//trim(char_field)//' line pattern 1'
    write (unit, *) '@    '//trim(char_field)//' baseline type 0'
    write (unit, *) '@    '//trim(char_field)//' baseline off'
    write (unit, *) '@    '//trim(char_field)//' dropline off'
    write (unit, *) '@    '//trim(char_field)//' fill type 0'
    write (unit, *) '@    '//trim(char_field)//' fill rule 0'
    write (unit, *) '@    '//trim(char_field)//' fill color 1'
    write (unit, *) '@    '//trim(char_field)//' fill pattern 1'
    write (unit, *) '@    '//trim(char_field)//' avalue off'
    write (unit, *) '@    '//trim(char_field)//' avalue type 2'
    write (unit, *) '@    '//trim(char_field)//' avalue char size 1.000000'
    write (unit, *) '@    '//trim(char_field)//' avalue font 0'
    write (unit, *) '@    '//trim(char_field)//' avalue color 1'
    write (unit, *) '@    '//trim(char_field)//' avalue rot 0'
    write (unit, *) '@    '//trim(char_field)//' avalue format general'
    write (unit, *) '@    '//trim(char_field)//' avalue prec 3'
    write (unit, *) '@    '//trim(char_field)//' avalue prepend ""'
    write (unit, *) '@    '//trim(char_field)//' avalue append ""'
    write (unit, *) '@    '//trim(char_field)//' avalue offset 0.000000 , 0.000000'
    write (unit, *) '@    '//trim(char_field)//' errorbar on'
    write (unit, *) '@    '//trim(char_field)//' errorbar place both'
    write (unit, *) '@    '//trim(char_field)//' errorbar color 2'
    write (unit, *) '@    '//trim(char_field)//' errorbar pattern 1'
    write (unit, *) '@    '//trim(char_field)//' errorbar size 1.000000'
    write (unit, *) '@    '//trim(char_field)//' errorbar linewidth 1.0'
    write (unit, *) '@    '//trim(char_field)//' errorbar linestyle 1'
    write (unit, *) '@    '//trim(char_field)//' errorbar riser linewidth 1.0'
    write (unit, *) '@    '//trim(char_field)//' errorbar riser linestyle 1'
    write (unit, *) '@    '//trim(char_field)//' errorbar riser clip off'
    write (unit, *) '@    '//trim(char_field)//' errorbar riser clip length 0.100000'
    write (unit, *) '@    '//trim(char_field)//' legend "'//trim(legend)//'"'

    ! write(unit,*) trim(char_field)//" hidden false"
    ! write(unit,*) trim(char_field)//" type xy"
    ! write(unit,*) trim(char_field)//" symbol 9"
    ! write(unit,*) trim(char_field)//" symbol size 1.000000"
    ! write(unit,*) trim(char_field)//" symbol color "//trim(char_colour)
    ! write(unit,*) trim(char_field)//" symbol pattern 1"
    ! write(unit,*) trim(char_field)//" symbol fill color 4"
    ! write(unit,*) trim(char_field)//" symbol fill pattern 1"
    ! write(unit,*) trim(char_field)//" symbol linewidth 1.5"
    ! write(unit,*) trim(char_field)//" symbol linestyle 1"
    ! write(unit,*) trim(char_field)//" symbol char 65"
    ! write(unit,*) trim(char_field)//" symbol char font 0"
    ! write(unit,*) trim(char_field)//" symbol skip 0"
    ! write(unit,*) trim(char_field)//" line type 0"
    ! write(unit,*) trim(char_field)//" line linestyle 1"
    ! write(unit,*) trim(char_field)//" lyine linewidth 1.0"
!s0 legend  ""
!s0 avalue on
!s0 avalue type 4
!s0 avalue char size 0.80000
!s0 avalue font 4
!s0 avalue color 1
!s0 avalue rot 0
!s0 avalue format general
!s0 avalue prec 4
!s0 avalue prepend ""
!s0 avalue append ""
!s0 avalue offset 0.000000 , 0.01000
  end subroutine xmgu_data_header

  !==================================================================!
  subroutine xmgu_vertical_line(unit, x_coord, y_max, y_min)
    !==================================================================!
    use od_constants, only: dp
    implicit none
    integer, intent(in)      :: unit

    real(dp), intent(in) :: x_coord, y_max, y_min

    character(20) ::  x_coord_char, y_max_char, y_min_char

    write (x_coord_char, '(F20.6)') x_coord
    write (y_max_char, '(F20.6)') y_max
    write (y_min_char, '(F20.6)') y_min

    x_coord_char = trim(adjustl((x_coord_char)))
    y_max_char = trim(adjustl((y_max_char)))
    y_min_char = trim(adjustl((y_min_char)))

    write (unit, *) '@with line'
    write (unit, *) '@    line on'
    write (unit, *) '@    line loctype world'
    write (unit, *) '@    line g0'
    write (unit, *) '@    line '//trim(x_coord_char)//', '//trim(y_min_char)&
     &//', '//trim(x_coord_char)//', '//trim(y_max_char)
    write (unit, *) '@    line linewidth 1.5'
    write (unit, *) '@    line linestyle 3'
    write (unit, *) '@    line color 2'
    write (unit, *) '@    line arrow 0'
    write (unit, *) '@    line arrow type 0'
    write (unit, *) '@    line arrow length 1.000000'
    write (unit, *) '@    line arrow layout 1.000000, 1.000000'
    write (unit, *) '@line def'

  end subroutine xmgu_vertical_line

end module xmgrace_utils
