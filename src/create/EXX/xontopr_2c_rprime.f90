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

 
! xontop_2c_rprime.f
! Program Description
! ===========================================================================
!       This routine calculates the integral over rprime, which is
! specifically used in the ontop interactions for the exact exchange.
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
        subroutine xontopr_2c_rprime (nspec_max, nssh, nalpha, itype,       &
    &                                 rcutoff, d, nrho, nz, lmax)
        use precision, only: wp
        use x_exact
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype
        integer, intent (in) :: lmax
        integer, intent (in) :: nalpha
        integer, intent (in) :: nrho
        integer, intent (in) :: nspec_max
        integer, intent (in) :: nz
 
        integer, intent (in) :: nssh (nspec_max)
 
        real(kind=wp), intent (in) :: d
        real(kind=wp), intent (in) :: rcutoff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer issh
        integer izp
        integer lqn
        integer nnz
 
        real(kind=wp) dzp
        real(kind=wp) psi3
        real(kind=wp) psi4
        real(kind=wp) r
        real(kind=wp) rho
        real(kind=wp) rhomax
        real(kind=wp) rhomin
        real(kind=wp) rp1
        real(kind=wp) rp2
        real(kind=wp) sumrp
        real(kind=wp) zmax
        real(kind=wp) zmin
        real(kind=wp) zp1
        real(kind=wp) zp2
 
        real(kind=wp), dimension (:), allocatable :: rhopmult
        real(kind=wp), dimension (:), allocatable :: zpmult
 
        real(kind=wp), external :: psiofr
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Set integration limits
        zmin = max(-rcutoff, d - rcutoff)
        zmax = min(rcutoff, d + rcutoff)
 
        rhomin = 0.0d0
        rhomax = rcutoff
 
! Strictly define what the density of the mesh should be. Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        dzp = rcutoff/real(nz)
        nnz = int((zmax - zmin)/dzp)
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1
 
        drhop = rcutoff/real(nrho)
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
 
! Loop over all of the shells of itype.
        allocate (rprime(nnrhop, 0:2*lmax, nssh(itype)))
        do issh = 1, nssh(itype)
 
! Loop over all the possible quantum numbers.
         do lqn = 0, 2*lmax
 
! Perform the radial integration over r' for each given r.
          do irho = 1, nnrhop
           r = rpoint(irho)
           if (r .lt. 1.0d-4) r = 1.0d-4
 
! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis). 
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
 
             if (rp2 .lt. rcutoff) then
              psi3 = psiofr (itype, nalpha, rp2)
              psi4 = psiofr (itype, issh, rp2)
 
! Limits from 0 to r.
              if (rp1 .le. r) then
               sumrp =                                                       &
      &         sumrp + zpmult(izp)*rhopmult(irhop)*rho*psi3*psi4            &
      &                            *rp1**lqn/r**(lqn + 1)
 
! Limits from r to rcutoff
              else
               sumrp =                                                       &
      &         sumrp + zpmult(izp)*rhopmult(irhop)*rho*psi3*psi4            &
      &                            *r**lqn/rp1**(lqn + 1)
              end if
             end if
            end do
           end do
 
! The actual integral as a function of r.
           rprime(irho, lqn, issh) = sumrp
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