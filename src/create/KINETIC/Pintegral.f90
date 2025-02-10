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

 
! Pintegral.f
! Program Description
! ===========================================================================
!
! This subroutine calculates the integral of two spherical harmonics of
! l,m and lp,mp multiplied by Legendre Polynomials of order from
! 0 to kmax.
!   integral (Ylm*Pn(cos(theta))Ylpmp)
!
! This will handle polynomials up to kmax = 10, in accordance
! with the subroutine integral which calculates the integral of two
! spherical harmonics multiplied by cos(theta) to the kth power.
!
! This will use the recursion relation for the Legendre Polynomials
!  (l+1)Pl+1(x) = (2l+1)*x*Pl(x) - l*Pl-1(x)
!   where x = cos(theta)
!
!  input:
!     l,m = l and m value for first spherical harmonic
!     lp,mp = l and m value for second spherical harmonic
!     kmax = highest order Legendre Polynomial to calculate
!           will calculate polynomials from 0 to kmax
!     Pint(i) = array of values for the integral of spherical harmonics
!          times Leg. Polynomial of order i
!
! ===========================================================================
! Code written by:
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
        subroutine Pintegral (l, m, lp, mp, kmax, Pint)
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer l
        integer m
        integer lp
        integer mp
        integer kmax
 
! Output
        real(kind=long) Pint (0:kmax)
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer nmax, imax, jmax
        parameter (nmax = 10, imax = 10, jmax=10)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i, j
 
        real(kind=long) arg
 
        real(kind=long) pn (0:imax,0:jmax)
        real(kind=long) theta (0:nmax)
 
        complex(kind=long) sum
 
! Procedure
! ===========================================================================
 
!----------------------------------------------------------------------
! Get the cosine integrals
!----------------------------------------------------------------------
        do i = 0, kmax
         call Tintegral (l, m, lp, mp, i, sum)
         theta(i) = dble(sum)
        end do
 
 
!----------------------------------------------------------------------
! Get the Legendre Polynomials by the recurrence relation
!----------------------------------------------------------------------
        do i = 0, imax
         do j = 0, jmax
          pn(i,j) = 0.0d0
         end do
        end do
 
        pn(0,0) = 1.0d0
        pn(1,1) = 1.0d0
 
        do i = 2, kmax
         do j = 0,i
          arg = dble(i)
          pn(i,j+1) = pn(i,j+1) + (2.0d0*arg - 1.0d0)*pn(i-1,j)/arg
          pn(i,j) = pn(i,j) - (arg - 1.0d0)*pn(i-2,j)/arg
         end do
        end do
 
 
!---------------------------------------------------------------------
! Calculate the integral of Legendre Polynomial with spherical
! harmonics.
!---------------------------------------------------------------------
        do i = 0, kmax
         Pint(i) = 0.0d0
         do j = 0, i
          Pint(i) = Pint(i) + pn(i,j)*theta(j)
         end do
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
