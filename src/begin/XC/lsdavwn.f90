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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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
        subroutine lsdavwn (rho, ex, ec, xpot, cpot)
        use constants
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real(kind=long), intent (in), dimension (2) :: rho

! Output
        real(kind=long), intent (out) :: ex
        real(kind=long), intent (out) :: ec

        real(kind=long), intent (out), dimension (2) :: cpot
        real(kind=long), intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: Ap = 0.0621814d0
        real(kind=long), parameter :: bp = 3.72744d0
        real(kind=long), parameter :: cp = 12.9352d0
        real(kind=long), parameter :: x0p = -0.10498d0

        real(kind=long), parameter :: Aa = 0.033773728d0
        real(kind=long), parameter :: ba = 1.13107d0
        real(kind=long), parameter :: ca = 13.0045d0
        real(kind=long), parameter :: x0a = -0.00475840

        real(kind=long), parameter :: Af = 0.01554535d0
        real(kind=long), parameter :: bf = 7.06042d0
        real(kind=long), parameter :: cf = 18.0578d0
        real(kind=long), parameter :: x0f = -0.32500d0

        real(kind=long), parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        real(kind=long) density, densitys

! spin polarization exchange-correlation pieces and derivatives
        real(kind=long) cte
        real(kind=long) d1ec, d2ec
        real(kind=long) d1ex, d2ex
        real(kind=long) ecA, ecAp
        real(kind=long) ecP, ecPp
        real(kind=long) ecF, ecFp
        real(kind=long) ecpx, ecpz
        real(kind=long) exP, exPp
        real(kind=long) expd, expz
        real(kind=long) jp, jf, ja
        real(kind=long) g, gp
        real(kind=long) h, hp
        real(kind=long) Qp, Qf, Qa
        real(kind=long) x, xp
        real(kind=long) XXp, XXf, XXa
        real(kind=long) zeta, zp1, zp2

! Allocate Arrays
! ===========================================================================

! Procedure
! =========================================================================
! Initialize some parameters
        density = rho(1) + rho(2)
        densitys = rho(1) - rho(2)
        zeta = densitys/density
        if (density .le. eps) then
         zeta = 0.0d0
         ec = 0.0d0
         ex = 0.0d0
         cpot(1) = 0.0d0
         cpot(2) = 0.0d0
         xpot(1) = 0.0d0
         xpot(2) = 0.0d0
         return
        end if

! Define simple derivatives
! *************************************************************************
        zp1 = 2.0d0*rho(1)/density**2
        zp2 = -2.0d0*rho(2)/density**2

        x = (3.0d0/(4.0d0*pi*density))**(1.0d0/6.0d0)
        xp = - (1.0d0/6.0d0)*x/density

        g = (9.0d0/8.0d0)*((1.0d0 + zeta)**(4.0d0/3.0d0)                     &
     &      + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0)
        gp = (3.0d0/2.0d0)*((1.0d0 + zeta)**(1.0d0/3.0d0)                    &
     &      - (1.0d0 - zeta)**(1.0d0/3.0d0))

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

        cte = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0) - 1.0d0))

        h = cte*((ecF - ecP)/ecA) - 1.d0
        hp = cte*((ecFp - ecPp)/ecA - (ecF - ecP)*(ecAp/ecA))

! Correlation functional ( and partials to z and x ):
        ec = ecP

        ecpx = ecPp
        ecpz = ecA*gp

! Partial derivatives VWN exchange functional
        d1ec = xp*ecpx + zp1*ecpz
        d2ec = xp*ecpx + zp2*ecpz

! ****************************************************************************
!
!       VNN EXCHANGE FUNCTIONAL
!
! ****************************************************************************
        exP = -(3.0d0/2.0d0)*(3.0d0*density/pi)**(1.0d0/3.0d0)
        exPp = exP/(3.0d0*density)

        ex = (1.0d0 + 4.0d0*g/9.0d0)*exP
        expd = ex/(3.0d0*density)
        expz = exP*gp

        d1ex = expd + zp1*expz
        d2ex = expd + zp2*expz

! ****************************************************************************
! Functions in Rydberg units - divide by factor of 2 to get Hartree
        ex = 0.5d0*ex
        ec = 0.5d0*ec

        xpot(1) = 0.5d0*density*d1ex + ex
        xpot(2) = 0.5d0*density*d2ex + ex
        cpot(1) = 0.5d0*density*d1ec + ec
        cpot(2) = 0.5d0*density*d2ec + ec

! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================

        return
        end
