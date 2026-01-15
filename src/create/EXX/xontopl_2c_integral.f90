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

 
! xontopl_2c_integral.f
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  The function V(1) is the
! non-local exchange potential, and is located at the site of one of the
! orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine xontopl_2c_integral (fraction, nalpha, lalpha, malpha,   &
     &                                  n1, l1, m1, n2, l2, m2, nz, nrho,   &
     &                                  d, itype1, itype2, rcutoff1,        &
     &                                  rcutoff2, lmax, sum)
        use precision, only: wp
        use coefficients
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype1
        integer, intent (in) :: itype2
 
        integer, intent (in) :: l1, l2      ! s, p, d, or f
        integer, intent (in) :: m1, m2      ! m quantum numbers
        integer, intent (in) :: n1, n2      ! which shell is being considered
        integer, intent (in) :: lalpha
        integer, intent (in) :: lmax
        integer, intent (in) :: malpha
        integer, intent (in) :: nalpha
        integer, intent (in) :: nrho        ! number of rho-points on grid
        integer, intent (in) :: nz          ! number of z-points on grid
 
        real(kind=wp), intent (in) :: d
        real(kind=wp), intent (in) :: fraction
        real(kind=wp), intent (in) :: rcutoff1
        real(kind=wp), intent (in) :: rcutoff2
 
! Output
        real(kind=wp), intent (out) :: sum
 
! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=wp), parameter :: eq2 = 14.39975d0
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer iz
        integer lqn
        integer mqn
        integer nnrho
        integer nnz
 
        real(kind=wp) cg1
        real(kind=wp) cg2
        real(kind=wp) drho
        real(kind=wp) dz
        real(kind=wp) factor
        real(kind=wp) phifactor
        real(kind=wp) psi1
        real(kind=wp) psi2
        real(kind=wp) r1, r2
        real(kind=wp) rho
        real(kind=wp) rhomax
        real(kind=wp) rhomin
        real(kind=wp) sumrp
        real(kind=wp) vofr
        real(kind=wp) z1, z2
        real(kind=wp) zmax
        real(kind=wp) zmin
 
        real(kind=wp), dimension (:), allocatable :: rhomult
        real(kind=wp), dimension (:), allocatable :: zmult
 
        real(kind=wp), external :: clebsch_gordon
        real(kind=wp), external :: delk
        real(kind=wp), external :: psiofr
        real(kind=wp), external :: rescaled_psi
        real(kind=wp), external :: rprimeofr
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Set integration limits
        zmin = max(-rcutoff1, d - rcutoff2)
        zmax = min(rcutoff1, d + rcutoff2)
 
        rhomin = 0.0d0
        rhomax = min(rcutoff1, rcutoff2)
 
! Strictly define what the density of the mesh should be. Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        dz = (rcutoff1 + rcutoff2)/(real(nz)*2.0d0)
        nnz = (zmax - zmin)/dz
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1
 
        drho = max(rcutoff1,rcutoff2)/real(nrho)
        nnrho = (rhomax - rhomin)/drho
        if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1
 
! Set up Simpson's rule factors. First for the z integration and then for
! the rho integration.
        allocate (zmult(nnz))
        zmult(1) = dz/3.0d0
        zmult(nnz) = dz/3.0d0
        do iz = 2, nnz - 1, 2
         zmult(iz) = 4.0d0*dz/3.0d0
        end do
        do iz = 3, nnz - 2, 2
         zmult(iz) = 2.0d0*dz/3.0d0
        end do
 
        allocate (rhomult(nnrho))
        rhomult(1) = drho/3.0d0
        rhomult(nnrho) = drho/3.0d0
        do irho = 2, nnrho - 1, 2
         rhomult(irho) = 4.0d0*drho/3.0d0
        end do
        do irho = 3, nnrho - 2, 2
         rhomult(irho) = 2.0d0*drho/3.0d0
        end do
 
! Initialize the sum to zero
        sum = 0.0d0
 
! ***************************************************************************
! The two-center integral calculated here is for a given matrix element, so
! a given lmu, mmu, and lnu, mnu.  Each matrix element involves a sum over
! all of the states at a given site, for the ontop (left) case this is a sum
! over all of the orbitals at the nu site.  The quantum number's which
! corresponding to one of these orbitals is passed as lalpha, malpha.
! Each lalpha, malpha combination are stored in a separate data file.
 
! Perform a sum over all quantum numbers l (up to lmax) and all
! corresponding quantum numbers m.
        do lqn = 0, 2*lmax
         do mqn = -lqn, lqn
 
! This factor comes from a result of multiplying the three Ylm's together.
! and after the factor of pi is multiplied out after the phi integration.
! Also, included in this factor is the Clebsch_Gordon coefficients which
! results from the angular integration over theta', phi'.
          cg1 = clebsch_gordon (lalpha, 0, lqn, 0, l1, 0)
          cg2 = clebsch_gordon (lalpha, malpha, lqn, mqn, l1, m1)
          phifactor = cg1*cg2*clm(l2,m2)*clm(lalpha,malpha)*clm(lqn,mqn)     &
     &                   *delk(m2,mqn+malpha)*sqrt(2.0d0*real(lalpha)+1.0d0) &
     &                   /(sqrt((2.0d0*real(lqn) + 1.0d0)                    &
     &                         *(2.0d0*real(l1) + 1.0d0))*2.0d0)
 
          if (phifactor .gt. 1.0d-4) then
 
! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis).
           do iz = 1, nnz
            z1 = zmin + real(iz-1)*dz
            z2 = z1 - d
            do irho = 1, nnrho
             rho = rhomin + real(irho-1)*drho
             r1 = sqrt(z1**2 + rho**2)
             r2 = sqrt(z2**2 + rho**2)
 
! Precaution against divide by zero
             if (r1 .lt. 1.0d-4) r1 = 1.0d-4
             if (r2 .lt. 1.0d-4) r2 = 1.0d-4
 
! Total integration factor
             factor = zmult(iz)*rhomult(irho)*phifactor
 
! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.
             if ((r1 .lt. rcutoff1) .and. (r2 .lt. rcutoff2)) then
              psi1 = psiofr (itype1, nalpha, r1)
              psi2 = psiofr (itype2, n2, r2)
 
! Perform the radial integration over r'.
! This integral was actually performed previous to calling this routine, since
! the integral answer can be performed independent of the r integration.
! Obtain the r dependent answer by interpolation.
              sumrp = rprimeofr (r1, lqn, n1, rcutoff1)
 
! Add magic factors based on what type of orbital is involved in the integration
! For the short-range coulomb interactions make spherically symmetric
              vofr = rescaled_psi (lqn, mqn, rho, r1, z1, 1.0d0)
              psi1 = rescaled_psi (lalpha, malpha, rho, r1, z1, psi1)
              psi2 = rescaled_psi (l2, m2, rho, r2, z2, psi2)
 
! This is the actual integral
              sum = sum + fraction*(eq2/2.0d0)*factor*psi1*vofr*psi2*sumrp*rho
             end if
 
! End loop for integration of z and rho.
            end do
           end do
          end if
 
! End loop over lqn and mqn
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
        deallocate (rhomult)
        deallocate (zmult)
 
! Format Statements
! ===========================================================================
 
        return
        end