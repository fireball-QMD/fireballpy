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


! dexc.f
! Program Description
! ===========================================================================
!
! This function computes
!
!    (n1+n2)*(exc(n1+n2) - vxc(n1+n2)) - sum_i ni*(exc(ni) - vxc(n1))
!
! with the densities ni of atom in_i at the distance r_i from their centers.
!
! On input:  in1,in2: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  dexc

! The catch comes in when we compute derivatives. We compute neutral,
! neutral for ideriv1. For other ideriv's we have the following KEY:
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
! ===========================================================================
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real*8 function dexc (in1, in2, r, z, r1, r2, iexc, fraction,
     1                        ix)
        implicit none

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1, in2
        integer ix

        real*8 fraction
        real*8 r
        real*8 r1
        real*8 r2
        real*8 z

! Local Parameters and Data Declaration
! ===========================================================================
        real*8 abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================
        integer in3

        real*8 dens
        real*8 densin
        real*8 densp
        real*8 denspin
        real*8 denspp
        real*8 densppin
        real*8 densz
        real*8 denszin
        real*8 denszz
        real*8 denszzin
        real*8 denspz
        real*8 denspzin
        real*8 dexc1c
        real*8 dexc2c
        real*8 dnuxc
        real*8 dnuxcs
        real*8 dnuxc2c
        real*8 dnuxcs2c
        real*8 dvxc1c
        real*8 dvxc2c
        real*8 hartree
        real*8 rin
        real*8 dexcc

! Procedure
! ===========================================================================
! By default set in3 = in2
        in3 = in2

! Two-center piece: [n1 + n2(r,z)]*(exc[n1 + n2(r,z)] - vxc[n1 + n2(r,z)])
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                      nnrho, nrho_points, nnz, nz_points, rho2c,
     2                      dens)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6
     1      .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points, rhop2c,
     2                       densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhozz2c, denszz)
        end if

! Convert to atomic units
        rin = r/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        denszin = densz*abohr**4
        densppin = denspp*abohr**5
        denspzin = denspz*abohr**5
        denszzin = denszz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc2c (iexc, fraction, rin, densin, denspin,
     1                    densppin, denszin, denszzin, denspzin, dexc2c,
     2                    dvxc2c, dnuxc2c, dnuxcs2c)

! Answers are in Hartrees convert to eV.
        hartree = 14.39975d0/abohr
        dexc = hartree*dens*(dexc2c - dvxc2c)

! One-center piece for atom 1: n1*(exc[n1(r1)] - vxc[n1(r1)])
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho,
     1                     dens, densp, denspp)

        rin = r1/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        densppin = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, densin, denspin,
     1           densppin, dexc1c, dvxc1c, dnuxc, dnuxcs, dexcc)

! Answers are in Hartrees convert to eV.
        dexc = dexc - hartree*dens*(dexc1c - dvxc1c)

! One-center piece for atom 2: n2*(exc[n2(r1)] - vxc[n2(r1)])
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in2
        call density_calc (iexc, ix, 2, in1, in2, in3, r2, drho,
     1                     dens, densp, denspp)

        rin = r2/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        densppin = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, densin, denspin,
     1          densppin, dexc1c, dvxc1c, dnuxc, dnuxcs, dexcc)

! Answers are in Hartrees convert to eV.
        dexc = dexc - hartree*dens*(dexc1c - dvxc1c)

! Format Statements
! ===========================================================================

        return
        end
