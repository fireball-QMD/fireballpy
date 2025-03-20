! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!            Jun Wang
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

!
! fireball-qmd is a free (GPLv3) open project.

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


! broadcast_MPI.f
! Program Description
! ===========================================================================
!       This file broadcasts the necessary information across the MPI
! processors. 
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah 
! 315 S. 1400 E. RM Dock
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353  
! Office telephone 801-585-1078 
! ===========================================================================
! 
! Program Declaration
! ===========================================================================
      	subroutine broadcast_MPI (signature)
      	implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input 
        character (len=70) signature

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
!  MPI Broadcast signature
        call MPI_BCAST (signature, 70, MPI_CHARACTER, 0, MPI_COMM_WORLD,     &
     &                  ierr_mpi)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

	return
        end
