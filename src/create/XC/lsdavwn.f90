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


! lsdavwn.f90
! Program Description
! ===========================================================================
!       This routine computes the exchange and correlation potenials and
! energies for the Vosko, Wilk, Nusair LSDA functional. Each spin component
! is considered.
!
! See
!      S.H. VOSKO and L. WILK and M. NUSAIR
!      Can. J. Phys., 58, 1200 (1980)
!
! ===========================================================================
! Code written by:
! Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine lsdavwn (rho, ex, ec, xpot, cpot, dnuxc, dnuxcs)
        use constants
        use precision, only: wp
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real(kind=wp), intent (in), dimension (2) :: rho

! Output
        real(kind=wp), intent (out) :: dnuxc
        real(kind=wp), intent (out) :: dnuxcs
        real(kind=wp), intent (out) :: ec
        real(kind=wp), intent (out) :: ex

        real(kind=wp), intent (out), dimension (2) :: cpot
        real(kind=wp), intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=wp), parameter :: Ap = 0.0621814d0
        real(kind=wp), parameter :: bp = 3.72744d0
        real(kind=wp), parameter :: cp = 12.9352d0
        real(kind=wp), parameter :: x0p = -0.10498d0

        real(kind=wp), parameter :: Aa = 0.033773728d0
        real(kind=wp), parameter :: ba = 1.13107d0
        real(kind=wp), parameter :: ca = 13.0045d0
        real(kind=wp), parameter :: x0a = -0.00475840d0

        real(kind=wp), parameter :: Af = 0.01554535d0
        real(kind=wp), parameter :: bf = 7.06042d0
        real(kind=wp), parameter :: cf = 18.0578d0
        real(kind=wp), parameter :: x0f = -0.32500d0

        real(kind=wp), parameter :: epsilon = 1.0d-10

! Local Variable Declaration and Description
! ===========================================================================
        real(kind=wp) density
        real(kind=wp) densitys

        real(kind=wp), dimension (3) :: cdpot
        real(kind=wp), dimension (3) :: xdpot

! spin polarization and derivatives
        real(kind=wp) zeta, zp1, zp2, zp1p2, zpp1, zpp2
        real(kind=wp) x, xp, xpp
        real(kind=wp) g, gp, gpp
        real(kind=wp) XXp, XXf, XXa , Qp, Qf, Qa, jp, jf, ja
        real(kind=wp) ecP, ecF, ecA, ecPp, ecFp, ecAp, ecPpp, ecFpp, ecApp
        real(kind=wp) cte, h, hp, hpp
        real(kind=wp) ecpx, ecpz, ecppx, ecppz, ecpxpz
        real(kind=wp) d1ec, d2ec, dd1ec, dd2ec, d1d2ec
        real(kind=wp) exP, exPp, exPpp
        real(kind=wp) expd, expz, exppd, exppz, expdpz
        real(kind=wp) d1ex, d2ex, dd1ex, dd2ex, d1d2ex
        real(kind=wp) F, Fs

! Allocate Arrays
! ===========================================================================

! Procedure
! =========================================================================
! Initialize some parameters
        stop ! zpp1 and zpp2 are not set
        density = rho(1) + rho(2)
        densitys = rho(1) - rho(2)
        zeta = densitys/density
        if (density .le. epsilon) then
         zeta = 0.0d0
         ec = 0.0d0
         ex = 0.0d0
         cpot = 0.0d0
         xpot = 0.0d0
         cdpot = 0.0d0
         xdpot = 0.0d0
         return
        end if

! Define simple derivatives
! *************************************************************************
        zp1 = 2.0d0*rho(1)/density**2
        zp2 = -2.0d0*rho(2)/density**2
        zp1p2 = 2.0d0*zeta/density**2

        x = (3.0d0/(4.0d0*pi*density))**(1.0d0/6.0d0)
        xp = - (1.0d0/6.0d0)*x/density
        xpp = 7.0d0*x/(36.0d0*density**2)

        g = (9.0d0/8.0d0)*((1.0d0 + zeta)**(4.0d0/3.0d0)                     &
     &     + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0)
        gp = (3.0d0/2.0d0)*((1.0d0 + zeta)**(1.0d0/3.0d0)                    &
     &      - (1.0d0 - zeta)**(1.0d0/3.0d0))
        gpp = (1.0d0/2.0d0)*((1.0d0 + zeta)**(-2.0d0/3.0d0)                  &
     &       - (1.0d0 - zeta)**(-2.0d0/3.0d0))

! Intermediate variables
        XXp = x**2.0d0 + bp*x + cp
        XXf = x**2.0d0 + bf*x + cf
        XXa = x**2.0d0 + ba*x + ca
        Qp = (4.0d0*cp - bp*bp)**0.5d0
        Qf = (4.0d0*cf - bf*bf)**0.5d0
        Qa = (4.0d0*ca - ba*ba)**0.5d0
        jp = 2.0d0*log(x - x0p) - log(XXp)                                   &
     &      + 2.0d0*((2.0d0*x0p + bp)/Qp)*atan(Qp/(2.0d0*x + bp))
        jf = 2.0d0*log(x - x0f) - log(XXf)                                   &
     &      + 2.0d0*((2.0d0*x0f + bf)/Qf)*atan(Qf/(2.0d0*x + bf))
        ja = 2.0d0*log(x - x0a) - log(XXa)                                   &
     &      + 2.0d0*((2.0d0*x0a + ba)/Qa)*atan(Qa/(2.0d0*x + ba))

! epsilon derivatives
        ecP = Ap*(2.0d0*log(x) - log(XXp)                                    &
     &            + (2.0d0*bp/Qp)*atan(Qp/(2.0d0*x + bp))                    &
     &            - (bp*x0p/(x0p*x0p + bp*x0p + cp))*jp)
        ecF = Af*(2.0d0*log(x) - log(XXf)                                    &
     &            + (2.0d0*bf/Qp)*atan(Qf/(2.0d0*x + bf))                    &
     &            - (bf*x0f/(x0f*x0f + bf*x0f + cf))*jp)
        ecA = Aa*(2.0d0*log(x) - log(XXa)                                    &
     &            + (2.0d0*ba/Qa)*atan(Qa/(2.0d0*x + ba))                    &
     &            - (ba*x0a/(x0a*x0a + ba*x0a + ca))*ja)

        ecPp = 2.0d0*Ap*cp/(XXp*x) - 2.0d0*Ap*bp*x0p/((x - x0p)*XXp)
        ecFp = 2.0d0*Af*cf/(XXf*x) - 2.0d0*Af*bf*x0f/((x - x0f)*XXf)
        ecAp = 2.0d0*Aa*ca/(XXa*x) - 2.0d0*Aa*ba*x0a/((x - x0a)*XXa)

        ecPpp = - 2.0d0*Ap*cp*(3.0d0*x**2 + 2.0d0*bp*x + cp)/(x*XXp)**2      &
     &          + 2.0d0*Ap*bp*x0p*((2.0d0*x + bp)*(x - x0p) + XXp)           &
     &                 /(XXp*(x - x0p))**2
        ecFpp = - 2.0d0*Af*cf*(3.0d0*x**2 + 2.0d0*bf*x + cf)/(x*XXf)**2      &
     &          + 2.0d0*Af*bf*x0f*((2.0d0*x + bf)*(x - x0f) + XXf)           &
     &                 /(XXf*(x - x0f))**2
        ecApp = - 2.0d0*Aa*ca*(3.0d0*x**2 + 2.0d0*ba*x + ca)/(x*XXa)**2      &
     &          + 2.0d0*Aa*ba*x0a*((2.0d0*x + ba)*(x - x0a) + XXa)           &
     &                 /(XXa*(x - x0a))**2

        cte = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0) - 1.0d0))

        h = cte*((ecF - ecP)/ecA) - 1.d0
        hp = cte*((ecFp - ecPp)/ecA - (ecF - ecP)*(ecAp/ecA))
        hpp = cte*((ecFpp - ecPpp)/ecA - (ecFp - ecPp)*ecAp/ecA**2           &
     &            - (ecF - ecP)*ecApp/ecA - (ecFp - ecPp)*ecAp/ecA           &
     &            + (ecF - ecP)*(ecAp/ecA)**2)


! Correlation functional ( and partials to z and x ):
        if (zeta .ne. 0.0d0) then
         ec = ecP + ecA*g*(1 + h*zeta**4)
        else
         ec = ecP
        end if

        ecpx = ecPp + ecAp*g*(1 + h*zeta**4) + ecA*g*hp*zeta**4
        ecpz = ecA*gp*(1.0d0 + h*zeta**4) + ecA*g*h*4*zeta**3

        ecppx = ecPp + ecApp*g*(1.0d0 + h*zeta**4) + 2.0d0*ecAp*g*hp*zeta**4 &
     &         + ecA*g*hpp*zeta**4
        ecppz = ecA*gpp*(1.0d0 + h*zeta**4) + ecA*gp*h*zeta**3               &
     &         + ecA*g*h*12.0d0*zeta**2
        ecpxpz = ecAp*gp*(1.0d0 + h*zeta**4) + ecA*gp*hp*zeta**4             &
     &          + ecAp*g*h*4.0d0*zeta**3 + ecA*g*hp*4.0d0*zeta**3

! Partial derivatives VWN exchanche functional
        d1ec = xp*ecpx + zp1*ecpz
        d2ec = xp*ecpx + zp2*ecpz

! Second partial derivatives
        dd1ec = xp**2*ecpp + 2.0d0*xp*zp1*ecpxpz + xpp*ecpx                  &
     &         + zp1*zp1*ecppz + zpp1*ecpz
        dd2ec = xp**2*ecpp + 2.0d0*xp*zp2*ecpxpz + xpp*ecpx                  &
     &         + zp2*zp1*ecppz + zpp2*ecpz
        d1d2ec = xp**2*ecpp+ xp*(zp1 + zp2)*ecpxpz + xpp*ecpx                &
     &          + zp1*zp2*ecppz + zp1p2*ecpz

! ****************************************************************************
!
!       VNN EXCHANGE FUNCTIONAL
!
! ****************************************************************************
        exP = (-3.0d0/2.0d0)*(3.0d0*density/pi)**(1.0d0/3.0d0)
        exPp = exP/(3.0d0*density)
        exPpp = - 2.0d0*exP/(3.0d0*density)**2

        ex =(1.0d0 + 4.0d0*g/9.0d0)*exP
        expd = ex/(3.0d0*density)
        exppd = -2.0d0*ex/(9.0d0*density**2)
        expz = exP*gp
        exppz = exP*gpp
        expdpz = exPp*gp

        d1ex = expd + zp1*expz
        d2ex = expd + zp2*expz

        dd1ex = exppd + 2.0d0*zp1*expdpz + expd + zp1*zp1*exppz + zpp1*expz
        dd2ex = exppd + 2.0d0*zp2*expdpz + expd + zp2*zp2*exppz + zpp2*expz
        d1d2ex = exppd + (zp1 + zp2)*expdpz + expd + zp1*zp2*exppz + zp1p2*expz

! Functions in Rydberg units - divide by factor of 2 to get Hartree
! ****************************************************************************
        xpot(1) = 0.5d0*(density*d1ex + ex)
        xpot(2) = 0.5d0*(density*d2ex + ex)
        cpot(1) = 0.5d0*(density*d1ec + ec)
        cpot(2) = 0.5d0*(density*d2ec + ec)
        ex = 0.5d0*ex
        ec = 0.5d0*ec

        cdpot(1) = 0.5d0*dd1ec
        cdpot(2) = 0.5d0*d1d2ec
        cdpot(3) = 0.5d0*dd2ec
        xdpot(1) = 0.5d0*dd1ex
        xdpot(2) = 0.5d0*d1d2ex
        xdpot(3) = 0.5d0*dd2ex

        dnuxc = 0.25d0*density*(xdpot(1) + 2.0d0*xdpot(2) + xdpot(3))        &
     &         + 0.5d0*(d1ec + d2ec) + 0.5d0*(d1ex + d2ex)                   &
     &         + 0.25d0*density*(cdpot(1) + 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.25d0*density*(xdpot(1) - 2.0d0*xdpot(2) + xdpot(3))       &
     &           + 0.5d0*(d1ec - d2ec) + 0.5d0*(d1ex - d2ex)                 &
     &           + 0.25d0*density*(cdpot(1) - 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.5d0*(ecA + 4.0d0*exP/9.0d0)/density


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end