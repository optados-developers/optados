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
!=========================================================================!
! MODULE od_pdis
! This module implements output of projected dispersion curves using
! routines from od_projection_utils.
!-------------------------------------------------------------------------!
module od_pdis

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  use od_constants, only: dp
  use od_projection_utils, only: projection_array, matrix_weights, max_am, proj_symbol, num_proj
  !-------------------------------------------------------------------------!

  private

  public :: pdis_calculate

contains

  subroutine pdis_calculate
    use od_electronic, only: elec_pdis_read, efermi_castep, efermi
    use od_projection_utils, only: projection_merge, projection_get_string, projection_analyse_orbitals
    use od_comms, only: on_root
    use od_parameters, only: set_efermi_zero, iprint
    use od_io, only: stdout

    implicit none

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                 Projected Dispersion Curve Calculation                     +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)')
    end if

    ! read in the pdos weights
    call elec_pdis_read

    ! look at the orbitals and figure out which atoms / states we have
    call projection_analyse_orbitals

    ! parse the pdis string to see what we want
    call projection_get_string

    ! form the right matrix elements
    call projection_merge

    ! set efermi
    if (on_root) write (stdout, '(1x,a1,a46,f8.4,a3,12x,a8)') "|",&
                          &" Set fermi energy from file : ", efermi_castep, " eV", "  <- EfC"
    efermi = efermi_castep

    ! write everything out
    if (on_root .and. (iprint > 2)) then
      call pdis_report_projectors
    end if

    if (on_root) then
      call pdis_write
    end if

  end subroutine pdis_calculate

  !===============================================================================
  subroutine pdis_write
    !===============================================================================
    ! Write out the pdis that was requested. Write them all to the same file, kpoint
    ! by kpoint.
    !===============================================================================
    use od_io, only: seedname
    implicit none
    character(len=512)  :: name

    ! write everything to one file
    name = trim(seedname)//'.pdis.dat'
    call write_pdis_to_file(1, num_proj, name)

  end subroutine pdis_write

  subroutine write_pdis_to_file(start_proj, stop_proj, name)
    !===============================================================================
    ! Write out projectors, start_proj, stop_proj, to file name, one kpoint at a time
    !===============================================================================
    use od_parameters, only: iprint, set_efermi_zero
    use od_algorithms, only: channel_to_am
    use od_electronic, only: pdos_mwab, all_kpoints, band_energy, efermi
    use od_cell, only: atoms_species_num, num_species, nkpoints
    use od_io, only: io_file_unit, io_error, io_date, stdout

    implicit none
    integer, intent(in) :: start_proj, stop_proj
    character(len=512), intent(in) :: name
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: string
    integer :: iproj, iam, ispecies_num, ispecies
    integer :: i, pdis_file, ierr
    integer :: N, n_eigen

    if (set_efermi_zero) then
      band_energy = band_energy - efermi
    end if

    write (string, '(I4,"(x,es14.7)")') (stop_proj - start_proj) + 1

    pdis_file = io_file_unit()
    open (unit=pdis_file, file=trim(name), iostat=ierr)
    if (iprint > 2) write (stdout, '(1x,a30,a30,17x,a1)') "| Writing PDIS projectors to: ", trim(name), "|"
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open output file in pdis: pdis_write")

    write (pdis_file, *) "##############################################################################"
    write (pdis_file, *) "#"
    write (pdis_file, *) "#                  O p t a D O S   o u t p u t   f i l e "
    write (pdis_file, '(1x,a1)') "#"
    call io_date(cdate, ctime)
    write (pdis_file, *) '#  Generated on ', cdate, ' at ', ctime
    write (pdis_file, '(1x,a78)') "##############################################################################"
    write (pdis_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
    write (pdis_file, '(1a,a)') '#', '|                    Projected Dispersion Curve -- Projectors                |'
    write (pdis_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'

    if (pdos_mwab%nspins > 1) then
      call io_error('pDIS not implemented for multiple spin channels.')
    else
      do iproj = start_proj, stop_proj
        write (pdis_file, '(1a,a1,a12,i4,a10,50x,a1)') '#', '|', ' Projector: ', iproj, ' contains:', '|'
        write (pdis_file, '(1a,a1,a16,10x,a14,36x,a1)') '#', '|', ' Atom ', ' AngM Channel ', '|'
        do ispecies = 1, num_species
          do ispecies_num = 1, atoms_species_num(ispecies)
            do iam = 1, max_am
              if (projection_array(ispecies, ispecies_num, iam, iproj) == 1) then
                write (pdis_file, '(1a,a1,a13,i3,a18,42x,a1)') "#", "|", proj_symbol(ispecies), &
                     &ispecies_num, channel_to_am(iam), '|'
              end if
            end do
          end do
        end do
        write (pdis_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
      end do

      do N = 1, nkpoints
        write (pdis_file, '(a10, i4, a10, es18.7, es18.7, es18.7)') 'K-point   ', N, '     ', (all_kpoints(i, N), i=1, 3)
        do n_eigen = 1, pdos_mwab%nbands
          write (pdis_file, '(es20.7,'//trim(string)//')') band_energy(n_eigen, 1, N), &
            (matrix_weights(i, n_eigen, N, 1), i=start_proj, stop_proj)
        end do
      end do

    end if

    close (pdis_file)

  end subroutine write_pdis_to_file

  subroutine pdis_report_projectors
    use od_algorithms, only: channel_to_am
    use od_cell, only: atoms_species_num, num_species
    use od_io, only: stdout
    implicit none

    integer :: iproj, iam, ispecies_num, ispecies

    write (stdout, *)
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    write (stdout, '(1x,a)') '|                    Projected Dispersion Curve -- Projectors                |'
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    do iproj = 1, num_proj
      write (stdout, '(1x,a1,a12,i4,a10,50x,a1)') '|', ' Projector: ', iproj, ' contains:', '|'
      write (stdout, '(1x,a1,a16,10x,a14,36x,a1)') '|', ' Atom ', ' AngM Channel ', '|'

      do ispecies = 1, num_species
        do ispecies_num = 1, atoms_species_num(ispecies)
          do iam = 1, max_am
            if (projection_array(ispecies, ispecies_num, iam, iproj) == 1) then
              write (stdout, '(1x,a1,a13,i3,a18,42x,a1)') "|", proj_symbol(ispecies), &
                ispecies_num, channel_to_am(iam), '|' !, " |  DEBUG :",  ispecies ,iam
            end if
          end do
        end do
      end do
      write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    end do
  end subroutine pdis_report_projectors

end module od_pdis
