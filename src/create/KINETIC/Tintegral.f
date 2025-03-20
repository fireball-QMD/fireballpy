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

 
! Tintegral.f
! Program Description
! ===========================================================================
! This subroutine will calculate the Tintegral of a product of two
! different spherical harmonics and cos(theta) to the power of k
!
!          integral ( Ylm*[cos(theta)]**kYl'm' )
!
! This routine uses this formula to expand the cos(theta)*Yl'm'
!
!      cos(theta)*Ylm = c1(l,m)*Yl+1,m + c2(l,m)*Yl-1,m
!
! where c1(l,m) = ((l-m+1)*(l+m+1)/(2*l+1)/(2*l+3))^1/2
!  and  c2(l,m) = ((l-m)(l+m)/(2l-1)/(2l+3))^1/2
!
! This is expanded using a tree algorithm as explained in the subroutine
!  tree.
!
! Once the (cos(theta))^k*Ylm is expanded, orthogonality relationships
!   are used to calculate the value of the integral.
!
! Right now the maximum k that can be done is k = 10
!   The number of array elements needed for any k is: 2**(k+1) - 1
!   so the cl and iyl arrays are dimensioned as cl(0:2046) and iyl(0:2046)
! =========================================================================
 
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
        subroutine Tintegral (l1, m1, l2, m2, k, sum)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer l1
        integer m1
        integer l2
        integer m2
        integer k
 
! Output:
        complex*16 sum             ! value of the integral
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer kmax
        integer maxnode
        parameter (kmax = 10, maxnode = 2046)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer nnodes
        integer istart
        integer iend
        integer i
        integer lnew
 
        integer iyl (0:maxnode)
 
        real*8 delk
 
        complex*16 cl (0:maxnode)
 
        external delk
 
! Procedure
! ===========================================================================
        if (k .gt. kmax) then
         write (*,*) ' The value of k',k,' given to routine integral '
         write (*,*) ' is too large, k must be less than or equal to '
     1               , kmax
         stop ' k > kmax in Tintegral.f '
        else if (k .lt. 0) then
         write (*,*) ' The value of k ', k, ' given to routine '
         write (*,*) ' integral is less than 0, k must be a positive '
         write (*,*) ' number for this subroutine. '
         stop ' k < 0 in Tintegral.f '
        end if
 
        nnodes = 2**(k + 1) - 1
        if (k .gt. 0) then
         call tree (l2, m2, nnodes, cl, iyl)
         sum = 0
         istart = 2**k - 1
         iend = 2*istart
         do i = istart, iend
          lnew = iyl(i)
          sum = sum + cl(i)*delk(l1,lnew)
         enddo
         sum = sum*delk(m1,m2)
        else
         sum = delk(l1,l2)*delk(m1,m2)
        endif
 
        return
        end
 
! Kronecker delta for l or m
! ===========================================================================
        function delk(i,j)
        implicit none
 
        integer i, j
        real*8 delk
 
        if (i .eq. j) then
         delk = 1.0d0
        else
         delk = 0.0d0
        end if
 
        return
        end
