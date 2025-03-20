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


! dpotxc12.f
! Program Description
! ===========================================================================
!
! This function computes the extended hubbard interaction: n1*n2*dnuxc(n1+n2)
! with the densities ni of atom in_i at the distance r_i from their centers.
!
! On input:  r, z: geometry information for the charge gradient
!
! On output:  dpotxc12
! ===========================================================================

! Code rewritten by:
! Otto F. Sankey (visiting)
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
        real*8 function dpotxc12s(r, z, iexc, fraction)
        implicit none

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real*8 fraction
        real*8 r
        real*8 z

! Local Parameters and Data Declaration
! ===========================================================================
        real*8 abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================

        real*8 dens
        real*8 densin
        real*8 densp
        real*8 denspp
        real*8 denspz
        real*8 densz
        real*8 denszz
        real*8 dexc2c
        real*8 dnuxc2c
        real*8 dnuxc2cs
        real*8 dvxc2c
        real*8 hartree
        real*8 rin

! Procedure
! ===========================================================================
! Two-center piece: [n1 + n2(r,z)]*(exc[n1 + n2(r,z)] - vxc[n1 + n2(r,z)])
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                      nnrho, nrho_points, nnz, nz_points, rho2c,
     2                      dens)

! Convert to atomic units
        rin = r/abohr
        densin = dens*abohr**3

! Here energy and potential due to exchange and correlation are calculated.
! OFS+JPL 1999 extended hubbard addition. We add nu to the output list.
! We only use this for interaction=12, and at this time only Ceperly Alder.
! Here we ONLY have Ceper. Alder. (See create.f)
        call get_potxc2c (iexc, fraction, rin, densin, densp, denspp,
     1                    densz, denszz, denspz, dexc2c, dvxc2c,
     2                    dnuxc2c, dnuxc2cs)

! dnuxc2c is energy*volume, so must be in Hartree*abohr**3
        hartree = 14.39975d0/abohr
        dpotxc12s = dnuxc2cs*hartree*abohr**3

! Format Statements
! ===========================================================================

        return
        end
