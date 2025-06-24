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
! MODULE od_pdos
! This module contains routines for writing the desired projector weights
! out to file, after accumulating across all kpoints.
!-------------------------------------------------------------------------------
module od_pdos

  !-------------------------------------------------------------------------!
  ! G L O B A L   V A R I A B L E S
  !-------------------------------------------------------------------------!
  use od_constants, only: dp
  use od_projection_utils, only: projection_array, matrix_weights, max_am, proj_symbol, num_proj, shortcut

  implicit none

  real(kind=dp), public, allocatable, save  :: dos_partial(:, :, :)

  !-------------------------------------------------------------------------!

  private

  public :: pdos_calculate

contains

  subroutine pdos_calculate
    use od_electronic, only: elec_pdos_read, efermi, efermi_set
    use od_dos_utils, only: dos_utils_calculate, dos_utils_set_efermi
    use od_projection_utils, only: projection_merge, projection_get_string, projection_analyse_orbitals
    use od_comms, only: on_root
    use od_parameters, only: iprint, set_efermi_zero
    use od_io, only: stdout

    implicit none

    if (on_root) then
      write (stdout, *)
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)') '+                 Projected Density Of States Calculation                    +'
      write (stdout, '(1x,a78)') '+============================================================================+'
      write (stdout, '(1x,a78)')
    end if

    ! read in the pdos weights
    call elec_pdos_read

    ! look at the orbitals and figure out which atoms / states we have
    call projection_analyse_orbitals

    ! parse the pdos string to see what we want
    call projection_get_string

    ! form the right matrix elements
    call projection_merge

    if (on_root .and. (iprint > 2)) then
      call pdos_report_projectors
    end if

    ! now compute the weighted dos
    call dos_utils_calculate(matrix_weights, dos_partial)

    ! and write everything out
    if (set_efermi_zero .and. .not. efermi_set) call dos_utils_set_efermi
    if (on_root) then
      call pdos_write
    end if

  end subroutine pdos_calculate

  !===============================================================================
  subroutine pdos_write
    !===============================================================================
    ! Write out the pdos that was requested. Write them all to the same file, unless
    ! we don't have a short cut. In this case, write 10 projectors per file.
    !===============================================================================
    use od_io, only: seedname
    implicit none

    character(len=20) :: start_iproj_name, end_iproj_name
    integer            ::  ifile, nfile, start_iproj, end_iproj
    character(len=512) :: name

    ! write(*,*) "======================================================================================"
    ! write(*,*) "ispecies,ispecies_num,iam,iproj,projection_array(ispecies,ispecies_num,iam,iproj)"
    ! do iproj=1,num_proj
    !    do iam=1,max_am
    !       do ispecies_num=1,maxval(atoms_species_num)
    !          do  ispecies=1,num_species
    !             write(*,*) ispecies,ispecies_num,iam,iproj,projection_array(ispecies,ispecies_num,iam,iproj)
    !          enddo
    !       enddo
    !    enddo
    ! enddo
    ! write(*,*) "======================================================================================"

    if (shortcut) then
      ! write everything to one file
      name = trim(seedname)//'.pdos.dat'
      call write_proj_to_file(1, num_proj, name)
    else ! not shortcut
      nfile = int(num_proj/10) + 1 ! Number of output files
      do ifile = 1, nfile
        start_iproj = (ifile - 1)*10 + 1 ! First projector in nfile
        if (ifile == nfile) then ! We're doing the last file
          end_iproj = num_proj
        else
          end_iproj = ifile*10
        end if

        write (start_iproj_name, '(I20.4)') start_iproj
        write (end_iproj_name, '(I20.4)') end_iproj
        name = trim(seedname)//'.pdos.proj-'//trim(adjustl(start_iproj_name))//'-'//trim(adjustl(end_iproj_name))//'.dat'
        call write_proj_to_file(start_iproj, end_iproj, name)
      end do
    end if
  end subroutine pdos_write

  subroutine write_proj_to_file(start_proj, stop_proj, name)
    !===============================================================================
    ! Write out projectors, start_proj, stop_proj, to file name
    !===============================================================================
    use od_dos_utils, only: E, dos_utils_set_efermi
    use od_parameters, only: dos_nbins, iprint, set_efermi_zero
    use od_algorithms, only: channel_to_am
    use od_electronic, only: pdos_mwab, efermi, efermi_set
    use od_cell, only: atoms_species_num, num_species
    use od_io, only: io_file_unit, io_error, io_date, stdout

    implicit none
    integer, intent(in) :: start_proj, stop_proj
    character(len=512), intent(in) :: name
    character(len=11) :: cdate
    character(len=9) :: ctime
    character(len=20) :: string
    integer :: iproj, iam, ispecies_num, ispecies
    integer :: idos, i, pdos_file, ierr
    real(kind=dp), allocatable :: E_shift(:)

    allocate (E_shift(dos_nbins), stat=ierr)
    if (ierr /= 0) call io_error('Error allocating E_shift in write_proj_to_file')
    if (set_efermi_zero) then
      E_shift = E - efermi
    else
      E_shift = E
    end if

    write (string, '(I4,"(1x,es14.7)")') (stop_proj - start_proj) + 1

    pdos_file = io_file_unit()
    open (unit=pdos_file, file=trim(name), iostat=ierr)
    if (iprint > 2) write (stdout, '(1x,a30,a30,17x,a1)') "| Writing PDOS projectors to: ", trim(name), "|"
    if (ierr .ne. 0) call io_error(" ERROR: Cannot open output file in pdos: pdos_write")

    write (pdos_file, '(a78)') "##############################################################################"
    write (pdos_file, '(1a)') "#"
    write (pdos_file, '(a57)') "#                  O p t a D O S   o u t p u t   f i l e "
    write (pdos_file, '(1x,a1)') "#"
    call io_date(cdate, ctime)
    write (pdos_file, '(a16, 11a, 4a, 9a)') '#  Generated on ', cdate, ' at ', ctime
    write (pdos_file, '(a78)') "##############################################################################"
    write (pdos_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
    write (pdos_file, '(1a,a)') '#', '|                    Partial Density of States -- Projectors                 |'
    write (pdos_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'

    if (pdos_mwab%nspins > 1) then
      do iproj = start_proj, stop_proj
        write (pdos_file, '(1a,a1,a12,i4,a10,50x,a1)') '#', '|', ' Column: ', iproj, ' contains:', '|'
        write (pdos_file, '(1a,a1,a16,10x,a14,5x,a15,16x,a1)') '#', '|', ' Atom ', ' AngM Channel ', ' Spin Channel ', '|'
        do ispecies = 1, num_species
          do ispecies_num = 1, atoms_species_num(ispecies)
            do iam = 1, max_am
            if (projection_array(ispecies, ispecies_num, iam, iproj) == 1) then
              write (pdos_file, '(1a,a1,a13,i3,a18,16x,a2,24x,1a)') "#", "|", proj_symbol(ispecies), &
                   &ispecies_num, channel_to_am(iam), 'Up', '|'
            end if
            end do
          end do
        end do
        write (pdos_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
      end do
      do iproj = start_proj, stop_proj
        write (pdos_file, '(1a,a1,a12,i4,a10,50x,a1)') '#', '|', ' Column: ', iproj + num_proj, ' contains:', '|'
        write (pdos_file, '(1a,a1,a16,10x,a14,5x,a15,16x,a1)') '#', '|', ' Atom ', ' AngM Channel ', ' Spin Channel ', '|'
        do ispecies = 1, num_species
          do ispecies_num = 1, atoms_species_num(ispecies)
            do iam = 1, max_am
              if (projection_array(ispecies, ispecies_num, iam, iproj) == 1) then
                write (pdos_file, '(1a,a1,a13,i3,a18,15x,a4,23x,1a)') "#", "|", proj_symbol(ispecies), &
                     &ispecies_num, channel_to_am(iam), 'Down', '|'
              end if
            end do
          end do
        end do
        write (pdos_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
      end do

      dos_partial(:, 2, :) = -dos_partial(:, 2, :)
      do idos = 1, dos_nbins
        write (pdos_file, '(es14.7,'//trim(string)//trim(string)//')') E_shift(idos), (dos_partial(idos, 1, i), &
             &i=start_proj, stop_proj), (dos_partial(idos, 2, i), i=start_proj, stop_proj)
      end do
    else
      do iproj = start_proj, stop_proj
        write (pdos_file, '(1a,a1,a12,i4,a10,50x,a1)') '#', '|', ' Column: ', iproj, ' contains:', '|'
        write (pdos_file, '(1a,a1,a16,10x,a14,36x,a1)') '#', '|', ' Atom ', ' AngM Channel ', '|'
        do ispecies = 1, num_species
          do ispecies_num = 1, atoms_species_num(ispecies)
            do iam = 1, max_am
              if (projection_array(ispecies, ispecies_num, iam, iproj) == 1) then
                write (pdos_file, '(1a,a1,a13,i3,a18,42x,a1)') "#", "|", proj_symbol(ispecies), &
                     &ispecies_num, channel_to_am(iam), '|'
              end if
            end do
          end do
        end do
        write (pdos_file, '(1a,a)') '#', '+----------------------------------------------------------------------------+'
      end do

      do idos = 1, dos_nbins
        write (pdos_file, '(es14.7,'//trim(string)//')') E_shift(idos), (dos_partial(idos, 1, i), i=start_proj, stop_proj)
      end do
    end if

    close (pdos_file)

    deallocate (E_shift, stat=ierr)
    if (ierr /= 0) call io_error('Error deallocating E_shift in write_proj_to_file')

  end subroutine write_proj_to_file

  subroutine pdos_report_projectors
    use od_algorithms, only: channel_to_am
    use od_cell, only: atoms_species_num, num_species
    use od_io, only: stdout
    implicit none

    integer :: iproj, iam, ispecies_num, ispecies

    write (stdout, *)
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    write (stdout, '(1x,a)') '|                    Partial Density of States -- Projectors                 |'
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    do iproj = 1, num_proj
      write (stdout, '(1x,a1,a12,i4,a10,50x,a1)') '|', ' Column: ', iproj, ' contains:', '|'
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
  end subroutine pdos_report_projectors

!!$!===============================================================================
!!$ subroutine count_atoms(orbital,num_orbitals,num_atoms)
!!$!===============================================================================
!!$! From the program LinDOS (AJM)
!!$! Take the orbial information and work out the number of atoms that the LCAO
!!$! describe
!!$!===============================================================================
!!$   use od_io, only : io_error
!!$   implicit none
!!$   integer, intent(in)            :: num_orbitals
!!$   type(orbitals), intent(inout)  :: orbital(1:num_orbitals) ! sepcies, ! num of each species ! l channel
!!$   integer, intent(out)           :: num_atoms
!!$
!!$   integer, allocatable           :: species_count(:)
!!$   integer                        :: num_species, ion_count
!!$   integer                        :: i, ierr
!!$
!!$   num_species=maxval(orbital(:)%species_no)  ! The maximum value is the highest species rank
!!$
!!$   allocate(species_count(1:num_species), stat=ierr)
!!$   if(ierr/=0) call io_error( " Error : cannot allocate species_count")
!!$
!!$   species_count=0
!!$   ion_count=0
!!$
!!$   do i=1,num_orbitals
!!$    ! If the species number is greater than the number we have for that species then use this
!!$    ! new number instead
!!$    ! NB I'm using data from the array orbital to index species count! :S
!!$    if(orbital(i)%rank_in_species>species_count(orbital(i)%species_no)) then
!!$      species_count(orbital(i)%species_no)=orbital(i)%rank_in_species
!!$      ion_count=ion_count+1
!!$    endif
!!$    orbital(i)%ion_no=ion_count
!!$   enddo
!!$
!!$   num_atoms=sum(species_count(:))
!!$
!!$   if(allocated(species_count)) then
!!$     deallocate(species_count, stat=ierr)
!!$     if(ierr/=0) stop " Error : cannot deallocate  species_count"
!!$   endif
!!$  end subroutine count_atoms

!!$  subroutine general_write_pdos
!!$    !===============================================================================
!!$    ! Write out the pdos that was requested. Make a pretty header so that the user
!!$    ! knows what each column means
!!$    !===============================================================================
!!$    use od_dos_utils,       only : E
!!$    use od_parameters,only : dos_nbins
!!$    use od_algorithms, only : channel_to_am
!!$    use od_electronic, only         : pdos_mwab
!!$    use od_cell, only : atoms_species_num, num_species
!!$    use od_io, only : io_file_unit, seedname, io_error, io_date
!!$
!!$    implicit none
!!$
!!$   character(len=11) :: cdate
!!$   character(len=9) :: ctime
!!$    character(len=20) :: string, filename
!!$    integer :: iproj, iam, ispecies_num, ispecies, species, species_num
!!$    integer :: last_species, last_species_num
!!$    integer :: idos, i, pdos_file,ierr, start_proj
!!$
!!$    logical :: projector_to_file
!!$
!!$
!!$    write(string,'(I4,"(x,es14.7)")') pdos_mwab%norbitals
!!$
!!$    start_proj=1
!!$    projectors: do iproj=1,num_proj
!!$       projector_to_file=.false.
!!$
!!$       ! Are we writing .pdos.projX.dat or .pdos.AtomAtomNo.dat?
!!$       ! does this projector contain more than one atom?
!!$       do iam=1,max_am
!!$          if(sum(projection_array(:,:,iam,iproj))>1) then
!!$             ! Yes it does contain more than one atom
!!$             projector_to_file=.true.
!!$          endif
!!$       enddo
!!$
!!$       if(projector_to_file) then
!!$          ! Then let's write out this projector and move on to the next one
!!$
!!$          ! Must first check whether this isn't the last one in a previous projector group
!!$          if(start_proj.ne.iproj) then ! Yes it is.
!!$             write(string,'(I20)') last_species_num
!!$             filename=proj_symbol(last_species)//adjustl(string)
!!$             write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ", &
!!$                  & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$             call write_proj_to_file(start_proj, iproj, filename)
!!$             start_proj=iproj+1 ! Reset start counter, and we've written the current one.
!!$             cycle projectors
!!$          endif
!!$
!!$          write(*,*) "For proj:", iproj, "there is more than one atom"
!!$          write(*,*) "Hence we're writing projectors to files"
!!$          write(string,'(I20)') iproj
!!$          filename="proj"//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ", &
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(iproj,iproj,filename)
!!$          start_proj=iproj+1 ! Reset start counter, and we've written the current one.
!!$          cycle projectors
!!$       endif
!!$
!!$       ! Since we're not writing just one projector to the file, we're going to have to work out
!!$       ! how many projectors there are going to be in this file.
!!$       ! Work out the species and species_rank of this projector
!!$       scan: do iam=1,max_am
!!$          do ispecies_num=1,maxval(atoms_species_num)
!!$             do  ispecies=1,num_species
!!$                if(projection_array(ispecies,ispecies_num,iam,iproj)==1) then
!!$                   species=ispecies
!!$                   species_num=ispecies_num
!!$                   write(*,*) "Projector ",iproj," is Species ", ispecies, " Rank ", ispecies_num
!!$                   exit scan
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo scan
!!$
!!$       ! First time through we just put the info about this projector into the registry
!!$       if(iproj==start_proj) then
!!$          write(*,*) "Skipping over projector:", iproj," as we've nothing to compare it against yet"
!!$          last_species=species
!!$          last_species_num=species_num
!!$          start_proj=iproj
!!$       ! If this is the same species as the last one. We go around again.
!!$       elseif((species==last_species).and.(species_num==last_species_num)) then
!!$          write(*,*) "Projector ", iproj, " has the same Species and Rank as ", start_proj
!!$          last_species=species
!!$          last_species_num=species_num
!!$       else ! We've come to the end of the projector group, so need to write out all the old ones.
!!$          write(*,*) "Projector ", iproj, " has different Species and Rank to ", start_proj
!!$          write(string,'(I20)') iproj
!!$          filename=trim(proj_symbol(species))//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj-1, " to file ", &
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(start_proj, iproj-1,filename)
!!$          start_proj=iproj ! Since we haven't written the current projector yet
!!$          last_species=species
!!$          last_species_num=species_num
!!$       endif
!!$
!!$       ! If this is our last loop, then we'd better write the last on out too.
!!$       if(iproj==num_proj) then
!!$          write(*,*) "Last Projector group ", start_proj, " to ", iproj
!!$          write(string,'(I20)') iproj
!!$          filename=trim(proj_symbol(species))//adjustl(string)
!!$          write(*,*) "So write out Projectors ", start_proj," to ", iproj, " to file ",&
!!$               & trim(seedname)//".pdos."//trim(filename)//".dat"
!!$          call write_proj_to_file(start_proj, iproj, filename)
!!$       endif
!!$
!!$    enddo projectors
!!$    return
!!$  end subroutine general_write_pdos

end module od_pdos
