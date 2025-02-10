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
!                      Richard B. Evans
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
! Texas A&M - Traian Dumitrica
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


! vxc.f
! Program Description
! ===========================================================================
!
! This function computes vxc(n1+n2) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  vxc: vxc[n1(r1) + n2(r2)]

! We calculate vxc(n1+n2). The catch comes in when we compute
! derivatives. We compute neutral, neutral for ideriv1. For other ideriv's we
! have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
!
! ======================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! James P. Lewis and Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ======================================================================
!
! Program Declaration
! ======================================================================
        real(kind=long) function vxc (r, z, iexc, fraction)
        use precision
        implicit none

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ======================================================================
! Input
        integer iexc

        real(kind=long) fraction
        real(kind=long) r
        real(kind=long) z

! Local Parameters and Data Declaration
! ======================================================================
        real(kind=long) abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ======================================================================
        real(kind=long) dens
        real(kind=long) densp
        real(kind=long) denspp
        real(kind=long) densz
        real(kind=long) denszz
        real(kind=long) denspz
        real(kind=long) dnuxc2c
        real(kind=long) dnuxcs2c
        real(kind=long) exc2c
        real(kind=long) hartree
        real(kind=long) rin
        real(kind=long) vxc2c

! Procedure
! ======================================================================
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                      nnrho, nrho_points, nnz, nz_points, rho2c, &
     &                      dens)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
     &      .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                       nnrho, nrho_points, nnz, nz_points, rhop2c, &
     &                       densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                       nnrho, nrho_points, nnz, nz_points, &
     &                       rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                       nnrho, nrho_points, nnz, nz_points, &
     &                       rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                       nnrho, nrho_points, nnz, nz_points, &
     &                       rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
     &                       nnrho, nrho_points, nnz, nz_points, &
     &                       rhozz2c, denszz)
        end if

! Convert to atomic units
        rin = r/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        densz = densz*abohr**4
        denspp = denspp*abohr**5
        denszz = denszz*abohr**5
        denspz = denspz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc2c (iexc, fraction, rin, dens, densp, denspp, &
     &                    densz, denszz, denspz, exc2c, vxc2c, dnuxc2c, &
     &                    dnuxcs2c)

! Answers are in Hartrees convert to eV.
        hartree = 14.39975d0/abohr
        vxc = hartree*vxc2c

! Format Statements
! ===========================================================================

        return
        end