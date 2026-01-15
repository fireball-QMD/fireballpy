! copyright info:
!
!                             @Copyright 1998
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
!                      Kirk VanOpdorp
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

 
! tree.f
! Program Description
! ===========================================================================
! This subroutine will expand the product of cos(theta)^k*Ylm
! where Ylm is a spherical harmonic for l and m.
! This is done by the formula:
!      cos(theta)*Yl,m = c1(l,m)*Yl+1,m + c2(l,m)*Yl-1,m
!  If one does this k times, a tree is made with each Yl,m getting two
!  offspring until you reach the desired k.
!
!                      k=0     Yl,m
!                             /    \
!             k=1  c1(l,m)Yl+1,m   c2(l,m)*Yl-1,m
!
! Each keeps branching doubly until the desired n.
!
! Note: As the branching continues, the coefficient of each child
!       is multiplied by the coefficient of the parent on down to
!       the end.
!
! Some useful formulas: The array starts at 0. The array index for the
! left child of indek i is 2*i+1. The array index for the right child
! of index i is 2*i+2.
! The array index for last parent: (n-2)/2.
! The array goes from 0 to n, where n is the number of nodes-1
!   nnodes = 2**(k+1) - 1
!
! The array index at the beginning of a row n, is given by:
!       2**n - 1
! The array index at the end of a row n, is given by:
!       2*(2**n - 1)
!
! input: l,m = l and m value for spherical harmonic being expanded
!          n = number of array elements.
! output: cl = array of values for coefficient
!        iyl = array of l-values for spherical harmonics being expanded
!
! Note: The arrays cl() and iyl() must be dimensioned in the calling
!       routine to be from zero to n or greater. Or we start overwriting
!       array values and get a bunch of garbage.
!
! ===========================================================================
! Code rewritten by:
! Kirk VanOpdorp
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-5909
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine tree (l, m, n, cl, iyl)
        use precision, only: wp
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer l
        integer m
        integer n               !number of nodes
 
! Output:
        integer iyl (0:n)       !array of l-values
 
        complex(kind=wp) cl(0:n)      !array of coefficients
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ifinal          !final parent
        integer i
        integer ileft           !left child
        integer iright          !right child
        integer lnow
 
        complex(kind=wp) c1           !function for c1
        complex(kind=wp) c2           !function for c2
 
        external c1, c2
 
! Procedure
! ===========================================================================
        cl(0) = 1
        iyl(0) = l
 
        ifinal = (n - 2)/2
 
        do i = 0, ifinal
         ileft = 2*i+1
         iright= 2*i+2
         lnow = iyl(i)
         cl(ileft) = c1(lnow,m)*cl(i)
         cl(iright) = c2(lnow,m)*cl(i)
         iyl(ileft) = iyl(i) + 1
         iyl(iright) = iyl(i) - 1
        end do
 
        return
        end
 
! Functions for calculating coefficients in angular integrals.
! ===========================================================================
        complex(kind=wp) function c1(l, m)
        use precision, only: wp
        implicit none
 
! Input
        integer l, m
 
        integer i, j
        complex(kind=wp) arg
 
        i = (l - m + 1)*(l + m + 1)
        j = (2*l + 1)*(2*l + 3)
        arg = cmplx(real(i, kind=wp)/real(j, kind=wp), kind=wp)
        c1 = sqrt(arg)
 
        return
        end
 
! ===========================================================================
        complex(kind=wp) function c2(l,m)
        use precision, only: wp
        implicit none
 
! Input
        integer l, m
 
        integer i, j
        complex(kind=wp) arg
 
        i = (l - m)*(l + m)
        j = (2*l - 1)*(2*l + 1)
        arg = cmplx(real(i, kind=wp)/real(j, kind=wp), kind=wp)
        c2 = sqrt(arg)
 
        return
        end