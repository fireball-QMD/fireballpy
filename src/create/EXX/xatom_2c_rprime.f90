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

 
! xatom_2c_rprime.f
! Program Description
! ===========================================================================
!       This routine calculates the integral over rprime, which is
! specifically used in the atom-atom interactions for the exact exchange.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine xatom_2c_rprime (n2, l2, m2, nalpha, lalpha, malpha,     &
    &                               itype1, itype2, rcutoff1, rcutoff2, d,  &
    &                               nrho, nz, lmax)
        use precision
        use x_exact
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype1
        integer, intent (in) :: itype2
        integer, intent (in) :: l2
        integer, intent (in) :: lalpha
        integer, intent (in) :: lmax
        integer, intent (in) :: m2
        integer, intent (in) :: malpha
        integer, intent (in) :: n2
        integer, intent (in) :: nalpha
        integer, intent (in) :: nrho
        integer, intent (in) :: nz
 
        real(kind=long), intent (in) :: d
        real(kind=long), intent (in) :: rcutoff1
        real(kind=long), intent (in) :: rcutoff2
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer izp
        integer lqn
        integer mqn
        integer nnz
 
        real(kind=long) dzp
        real(kind=long) psi3
        real(kind=long) psi4
        real(kind=long) r
        real(kind=long) rho
        real(kind=long) rhomax
        real(kind=long) rhomin
        real(kind=long) rp1
        real(kind=long) rp2
        real(kind=long) sumrp
        real(kind=long) vofr
        real(kind=long) zmax
        real(kind=long) zmin
        real(kind=long) zp1
        real(kind=long) zp2
 
        real(kind=long), dimension (:), allocatable :: rhopmult
        real(kind=long), dimension (:), allocatable :: zpmult
 
        real(kind=long), external :: psiofr
        real(kind=long), external :: rescaled_psi
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Set integration limits
        zmin = max(-rcutoff1, d - rcutoff2)
        zmax = min(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = min(rcutoff1,rcutoff2)

! Strictly define what the density of the mesh should be. Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        dzp = (rcutoff1 + rcutoff2)/real(2*nz)
        nnz = int((zmax - zmin)/dzp)
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1

        drhop = max(rcutoff1,rcutoff2)/real(nrho)
        nnrhop = int((rhomax - rhomin)/drhop)
        if (mod(nnrhop,2) .eq. 0) nnrhop = nnrhop + 1

! Set up Simpson's rule factors. First for the rho integration and then for
! the z integration.
        allocate (rpoint(nnrhop))
        allocate (rhopmult(nnrhop))
        rpoint(1) = rhomin
        rpoint(nnrhop) = rhomax
        rhopmult(1) = drhop/3.0d0
        rhopmult(nnrhop) = drhop/3.0d0
        do irho = 2, nnrhop - 1, 2
         rpoint(irho) = rhomin + real(irho - 1)*drhop
         rhopmult(irho) = 4.0d0*drhop/3.0d0
        end do
        do irho = 3, nnrhop - 2, 2
         rpoint(irho) = rhomin + real(irho - 1)*drhop
         rhopmult(irho) = 2.0d0*drhop/3.0d0
        end do
 
        allocate (zpmult(nnz))
        zpmult(1) = dzp/3.0d0
        zpmult(nnz) = dzp/3.0d0
        do izp = 2, nnz - 1, 2
         zpmult(izp) = 4.0d0*dzp/3.0d0
        end do
        do izp = 3, nnz - 2, 2
         zpmult(izp) = 2.0d0*dzp/3.0d0
        end do

! Loop over all the possible quantum numbers.
        allocate (rprime(nnrhop, 0:2*lmax, -2*lmax:2*lmax))
        do lqn = 0, 2*lmax
         do mqn = -lqn, lqn
 
! Perform the radial integration over r' for each given r.
          do irho = 1, nnrhop
           r = rpoint(irho)
           if (r .lt. 1.0d-4) r = 1.0d-4

! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis). We are integrating in cylindrical coordinates
! here so the phi integration only goes from 0 - pi, rather than 0 - 2*pi;
! therefore, we divide the rho by a factor of two here.  
           sumrp = 0.0d0
           do izp = 1, nnz 
            zp1 = zmin + real(izp - 1)*dzp
            zp2 = zp1 - d
            do irhop = 1, nnrhop
             rho = rpoint(irhop)
             rp1 = sqrt(zp1**2 + rho**2)
             rp2 = sqrt(zp2**2 + rho**2)

! Precaution against divide by zero
             if (rp1 .lt. 1.0d-4) rp1 = 1.0d-4
             if (rp2 .lt. 1.0d-4) rp2 = 1.0d-4

             if (rp1 .lt. rcutoff1 .and. rp2 .lt. rcutoff2) then
              psi3 = psiofr (itype2, nalpha, rp2)
              psi4 = psiofr (itype1, n2, rp1)
 
! Add magic factors based on what type of orbital is involved in the integration
! For the short-range coulomb interactions make spherically symmetric
              vofr = rescaled_psi (lqn, mqn, rho, rp1, zp1, 1.0d0)
              psi3 = rescaled_psi (lalpha, malpha, rho, rp2, zp2, psi3)
              psi4 = rescaled_psi (l2, m2, rho, rp1, zp1, psi4)

! Limits from 0 to r.
              if (rp1 .le. r) then
               sumrp = sumrp + zpmult(izp)*rhopmult(irhop)*rho*psi3*psi4*vofr &
      &                                   *rp1**lqn/r**(lqn + 1)
 
! Limits from r to rcutoff
              else
               sumrp = sumrp + zpmult(izp)*rhopmult(irhop)*rho*psi3*psi4*vofr &
      &                                   *r**lqn/rp1**(lqn + 1)
              end if
             end if
            end do
           end do
 
! The actual integral as a function of r.
           rprime(irho, lqn, mqn) = sumrp
          end do
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
        deallocate (rhopmult)
        deallocate (zpmult)
 
! Format Statements
! ===========================================================================
 
        return
        end