! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

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


! get_potxc.f90
! Program Description
! ===========================================================================
!  This program will access the requested exchange and correlation
!  functionals from the following list:
!          1  LDA   Wigner
!          2  LDA   Hedin/Lundqvist
!          3  LDA   Ceperley/Alder Perdew/Zunger (1980)
!          4  GGA   Perdew/Wang (1991)
!          5  GGA   Becke (1988) X, Perdew (1986) C
!          6  GGA   Perdew/Burke/Ernzerhof (1996)
!          7  LDA   Zhao/Parr
!          8  LDA   Ceperley/Alder Perdew/Wang (1991)
!          9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!         10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
!         11  LSDA  Volko/Wilk/Nusair (1980)
! The numerical value above is assigned to the variable iexc output
! exchange-correlation potential. This program has been modified to account
! for the fact that the density is a sum of two densities at two different
! centers.  Also, the potential is evaluated at one point in space and the
! integrals (matrix elements) are evaluated elsewhere.
!
! input
!    iexc        xc scheme
!    r           radial coordinate
!    rho         sum of atomic density in au
!    rhop        sum of atomic density gradient (with respect to r) in au
!    rhopp       sum of second gradient (with respect to r) in au
!    rhoz        sum of atomic density gradient (with respect to z) in au
!    rhozz       sum of second gradient (with respect to z) in au
!
! output
!    vpxc        xc potential
!    newexc      xc energy
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Eduardo Mendez
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
        subroutine get_potxc2c (iexc, fraction, r, rho, rhop, rhopp, rhoz,  &
     &                          rhozz, rhopz, newexc, vpxc, dnuxc, dnuxcs)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real(kind=long), intent(in) :: fraction
        real(kind=long), intent(inout) :: r
        real(kind=long), intent(inout) :: rho
        real(kind=long), intent(in) :: rhop
        real(kind=long), intent(in) :: rhopp
        real(kind=long), intent(in) :: rhopz
        real(kind=long), intent(in) :: rhoz
        real(kind=long), intent(in) :: rhozz

! Output
        real(kind=long), intent(out) :: dnuxc
        real(kind=long), intent(out) :: dnuxcs
        real(kind=long), intent(out) :: newexc
        real(kind=long), intent(out) :: vpxc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ix

        real(kind=long) aln
        real(kind=long) dec
        real(kind=long) dex
        real(kind=long) drvexc
        real(kind=long) ecp
        real(kind=long) ex
        real(kind=long) exc
        real(kind=long) fx
        real(kind=long) fxc
        real(kind=long) rs
        real(kind=long) x
        real(kind=long) zeta

        real(kind=long), dimension (2) :: cpot
        real(kind=long), dimension (2) :: d
        real(kind=long), dimension (2) :: dp
        real(kind=long), dimension (2) :: dpp
        real(kind=long), dimension (2) :: dpz
        real(kind=long), dimension (2) :: dz
        real(kind=long), dimension (2) :: dzz
        real(kind=long), dimension (2) :: xpot

! Procedure
! ===========================================================================
! Initialize variables (only used in 3,11)
        dnuxc = 0.0d0
        dnuxcs = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
         rho = 0.0d0
         dnuxc = 0.0d0
         newexc = 0.0d0
         vpxc = 0.0d0
         return
        else if (rho .lt. 1.0d-5) then
         rho = 1.0d-5
        end if

! Determine exchange-correlation potentials lsdavwn
        if (iexc .eq. 11) then
         zeta = 0.0d0
         d(1) = 0.5d0*rho*(1 + zeta)
         d(2) = 0.5d0*rho*(1 - zeta)
         call lsdavwn (d, dex, dec, xpot, cpot, dnuxc, dnuxcs)
         vpxc = xpot(1) + cpot(1)                ! Holds only for zeta = 0 only

! correlation (C) Wigner
        else if (iexc .eq. 1) then
         call wigner (rho, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! C Hedin-Lundqvist
        else if (iexc .eq. 2) then
         if (rho .ne. 0.0d0) then
          rs = 0.62035049d0*rho**(-1.0d0/3.0d0)
          x = rs/21.0d0
          aln = log(1.0d0 + 1.0d0/x)
          ecp = aln + (x**3*aln - x*x) + x/2 - 1.0d0/3.0d0
          dex = -0.458175d0/rs - 0.0225d0*ecp
          vpxc = -0.6109d0/rs - 0.0225d0*aln
         else
          dex = 0.0d0
          vpxc = 0.0d0
         end if

! XC Ceperley - Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
        else if (iexc .eq. 3) then
         call ceperley_alder (rho, ex, fx, exc, fxc, drvexc, dnuxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! X Perdew, C Perdew, generalized gradient approximation 1992
        else if (iexc .eq. 4) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (3, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (2, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! X Becke, C Perdew, generalized gradient approximation
        else if (iexc .eq. 5) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (2, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (3, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! XC Wigner-scaled LDA of PRA 46, R5320 (1992)
        else if(iexc .eq. 7) then
         call wigscaled (rho, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! XC Ceperley-Alder in Perdew-Wang parametrization of 1991
        else if (iexc .eq. 8) then
         d = 0.5d0*rho
         dp = 0.0d0
         dpp = 0.0d0
         dz = 0.0d0
         dzz = 0.0d0
         dpz = 0.0d0
         call ggaxrad2c (1, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (1, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! C Lee-Yang-Parr
        else if (iexc .eq. 9 .or. iexc .eq. 10 .or. iexc .eq. 12) then

! X Becke gga by default
         ix = 2

! X Perdew-Wang gga
         if (iexc .eq. 10) ix = 3
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (ix, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (4 , r, d, dp, dpp, dz, dzz,      cpot, dec)
         if (iexc .ne. 12) then
          vpxc = xpot(1) + cpot(1)
         else
          vpxc = (1.0d0 - fraction)*xpot(1) + cpot(1)
          dex = (1.0d0 - fraction)*dex
         end if

! XC burke-perdew-ernzerhof gga 1996
        else if(iexc .eq. 6) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (5, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (5, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
         write (*,*) ' In get_potxc2c.f90 - '
         write (*,*) ' stop: xc option not implemented', iexc
         stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dec + dex

! Format Statements
! ===========================================================================

        return
        end