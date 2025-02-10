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

 
! clebsch_gordon.f
! Program Description
! ===========================================================================
!       This routine calculates the Clebsch-Gordon coefficients which
! are represented by - <l1,l2;m1,m2|l,m>.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real(kind=long) function clebsch_gordon (l1, m1, l2, m2, l, m)
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: l
        integer, intent (in) :: l1
        integer, intent (in) :: l2
        integer, intent (in) :: m
        integer, intent (in) :: m1
        integer, intent (in) :: m2
 
! Local Parameters and Data Declaration
! ===========================================================================
! The maximum z value is to be only 6, since the lmax value is only 6 in the
! routine which calls this function.
        integer, parameter :: izmax = 6
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iz
 
        real(kind=long) piece1
        real(kind=long) piece2
        real(kind=long) piece3
 
        real(kind=long), external :: factorial
 
! Procedure
! ===========================================================================
! Initialize the coefficient to zero.
        clebsch_gordon = 0.0d0
 
! First determine that the values of the l1, l2, and l satisfiy the
! triangle equation - |l1-l2| =< l >= l1 + l2.
        if ((m .eq. 0) .and. (m1 .eq. 0) .and. (m2 .eq. 0) .and.           &
     &      (mod(l + l1 + l2, 2) .ne. 0)) return
        if (l .lt. abs(l1 - l2)) return
        if (l .gt. (l1 + l2)) return
 
! The other condition is that m = m1 + m2
        if (m .ne. (m1 + m2)) return
 
! The clebsch_gordon coefficient will be written as a product of three
! pieces. So clebsch_gordon = piece1*piece2*piece3
        piece1 = sqrt(((2.0d0*l + 1)*factorial(l1 + l2 - l)                  &
     &                 *factorial(l1 - l2 + l)*factorial(-l1 + l2 + l))      &
     &                /factorial(l1 + l2 + l + 1))
 
        piece2 = sqrt(factorial(l1 + m1)*factorial(l1 - m1)                  &
     &                *factorial(l2 + m2)*factorial(l2 - m2)                 &
     &                *factorial(l + m)*factorial(l - m))
 
        piece3 = 0.0d0
        do iz = 0, izmax
         if (((l1 + l2 - l - iz) .ge. 0) .and. ((l1 - m1 - iz) .ge. 0)       &
     &       .and. ((l2 + m2 - iz) .ge. 0) .and. ((l - l2 + m1 + iz) .ge. 0) &
     &       .and. ((l - l1 - m2 + iz) .ge. 0)) then
          piece3 = piece3 + (-1)**iz                                         &
     &                /(factorial(iz)*factorial(l1 + l2 - l - iz)            &
     &                *factorial(l1 - m1 - iz)*factorial(l2 + m2 - iz)       &
     &                *factorial(l - l2 + m1 + iz)*factorial(l - l1 - m2 + iz))
         end if
        end do
 
! Now calculate the coefficient
        clebsch_gordon = piece1*piece2*piece3
 
! Format Statements
! ===========================================================================
 
        return
        end