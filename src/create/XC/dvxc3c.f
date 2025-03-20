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


! dvxc3c.f
! Program Description
! ===========================================================================
!
!       This subroutine computes vxc(n1+n2+n3) - vxc(n1+n2) with the
! densities ni of atom in_i at the distance ri from their centers.  Due to
! the small contribution of the three center case to the overall energy,
! only the lda level of theory will be used for this calculation.
!
! Input Variables
!
!     iexc:        Desired exchange/correlation functional
!     in1,in2,in3: Atomic indices
!     r1,r2,r3:    Radii for n1, n2, and n3
!     IX:          Derivative w.r.t the j1-th charge using the
!                  j2-th magnitude
!
! Common Variables
!
!     orbocc:      Orbital occupations
!     numorb:      Number of orbitals
!
! Output Variables
!
!     dvxc3c:      vxc(n1+n2+n3) - vxc(n1+n2)
!
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! Richard B. Evans
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
        real*8 function dvxc3c (iexc, r1, r2, r3, in1, in2, in3, IX)
        implicit none

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1
        integer in2
        integer in3
        integer ix

        real*8 r1
        real*8 r2
        real*8 r3

! Local Parameters and Data Declaration
! ===========================================================================
        real*8 abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================

        real*8 dens
        real*8 dens1
        real*8 dens2
        real*8 dens3
        real*8 densp
        real*8 densp1
        real*8 densp2
        real*8 densp3
        real*8 denspp
        real*8 denspp1
        real*8 denspp2
        real*8 denspp3
        real*8 dnuxc
        real*8 dnuxcs
        real*8 exc2c
        real*8 exc3c
        real*8 fraction
        real*8 hartree
        real*8 rin
        real*8 vxc2c
        real*8 vxc3c
        real*8 dexc

! Procedure
! ===========================================================================
! This exchange-correlation routine deals only with three-center interactions,
! therefore, we do not do exact exchange for the three-center terms and
! fraction should be initiallized to 1.0d0.
        fraction = 1.0d0

! Evaluate the spherically averaged wave functions for all orbitals of the
! two atoms
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho,
     1                     dens1, densp1, denspp1)
        call density_calc (iexc, ix, 2, in1, in2, in3, r2, drho,
     1                     dens2, densp2, denspp2)
        call density_calc (iexc, ix, 3, in1, in2, in3, r3, drho,
     1                     dens3, densp3, denspp3)

! Three-center-piece: dvxc3c[n1(r1) + n2(r2) + n3(r3)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total three-center density is the sum of the three.
! Note that rin, densp, and denspp are not used in the LDA limits.
        dens = dens1 + dens2 + dens3
        densp = densp1 + densp2 + densp3
        denspp = denspp1 + denspp2 + denspp3

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp,
     1                    exc3c, vxc3c, dnuxc, dnuxcs, dexc)

! Two-center-piece: dvxc2c[n1(r1) + n2(r2)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total two-center density is the sum of the two - dens1 + dens2.
! Note that rin, densp, and denspp are not used in the LDA limits.
        dens = dens1 + dens2
        densp = densp1 + densp2
        denspp = denspp1 + denspp2

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp,
     1                    exc2c, vxc2c, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        hartree = 14.39975d0/abohr
        dvxc3c = hartree*(vxc3c - vxc2c)

! Format Statements
! ===========================================================================

        return
        end
