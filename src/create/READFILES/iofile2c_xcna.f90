! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! University of Utah - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold

!
! fireball-qmd is a free (GPLv3) open project.

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


! iofile2c_xcna.f90
! Program Description
! ===========================================================================
!       This subroutine takes a given root and suffix with three index 
! numbers.  This information is then combined to give an input or output 
! filename.  This file is then opened according to the input unit device 
! given.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine iofile2c_xcna (root, suffix, isorp, index1, index2, iunit,&
     &                            filename, skip)
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: index1
        integer, intent (in) :: index2
        integer, intent (in) :: isorp
        integer, intent (in) :: iunit
 
        character (len=20), intent (in) :: root
        character (len=3), intent (in) :: suffix
 
! Output
        character (len=40), intent (out) :: filename
 
        logical, intent (out) :: skip

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nplace = 2
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i, j
        integer lr
        integer ls
 
        character (len=80) form
 
! Procedure
! ===========================================================================
! Find true length of root and suffix
        lr = len(root)
        do i = lr, 1, -1
         if (root(i:i) .ne. ' ') exit 
        end do
        lr = i
        ls = len(suffix)
        do i = ls, 1, -1
         if (suffix(i:i) .ne. ' ') exit
        end do
        ls = i
        if (lr + ls + nplace + 1 .gt. len(filename)) then
         write (*,'('' filename too short for input in iofile'')')
         write (*,*) root, suffix, isorp, index1, index2, iunit
         stop ' error in iofile2c_xcna.f90 '
        end if
 
! Form file name
! ***************************************************************************
! Fill in root here
        do i = 1, lr
         filename(i:i) = root(i:i)
        end do
 
! Test to make sure index is not too big
        if (index1 .lt. 0 .or. index2 .lt. 0 .or. isorp .lt. 0 .or.          &
     &      index1 .ge. 10**nplace .or. index2 .ge. 10**nplace) then
         write (*,'('' index out of range in iofile'')')
         write (*,*) root, suffix, isorp, index1, index2, iunit
         stop ' error in iofile2c_xcna.f90 '
        end if
 
! write the format to form this will be (i3.3) for nplace = 3
        write (form,'(''(i'',i1,''.'',i1,'')'')') nplace, nplace
 
! Write out the index numbers with a '.' between each number
        filename(lr+1:lr+1) = '_'
        write (filename(lr+2:lr+nplace+1),form) isorp
        filename(lr+nplace+2:lr+nplace+2) = '.'
        write (filename(lr+nplace+3:lr+2*nplace+2),form) index1
        filename(lr+2*nplace+3:lr+2*nplace+3) = '.'
        write (filename(lr+2*nplace+4:lr+3*nplace+3),form) index2
        filename(lr+3*nplace+4:lr+3*nplace+4) = '.'
 
! Fill in the suffix here
        do i = 1, ls
         j = lr + 3*nplace + 4 + i
         filename(j:j) = suffix(i:i)
        end do
 
! Fill rest with spaces in case there is garbage there now
        do i = lr + 3*nplace + ls + 5, len(filename)
         filename(i:i) = ' '
        end do
 
! Filename is formed, open the file
        skip = .false.
        inquire (file = filename, exist = skip)
        open (unit = iunit, file = filename, status = 'unknown')

! Format Statements
! ===========================================================================
 
        return
        end