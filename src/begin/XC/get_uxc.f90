! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! get_uxc.f90
! Program Description
! ===========================================================================
!
!       This routine calculates the exchange-correlation potential and
! energies for a variety of different functional forms. The variable ioption
! is input to direct which functional form is to be used. The options thus
! far are:
!
!        1   LDA  Wigner
!        2   LDA  Hedin/Lundqvist
!        3   LDA  Ceperley/Alder Perdew/Zunger (1980)
!        4   GGA  Perdew/Wang (1991)
!        5   GGA  Becke (1988) X, Perdew (1986)
!        6   GGA  Perdew/Burke/Ernzerhof (1996)
!        7   LDA  Zhao/Parr
!        8   LDA  Ceperley/Alder Perdew/Wang (1991)
!        9   GGA  Becke (1988) X, Lee/Yang/Parr (1988)
!        10  GGA  Perdew/Wang (1991) X, Lee/Yang/Parr (1988)
!        11  LSDA Volko/Wilk/Nusair (1980)
!        12  B3LYP  mixing exact exchange and BLYP (9 GGA)
!
! The original code is from
! Martin Fuchs, FHI der MPG, Berlin, 01-1993
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_uxc (ioption, mesh, r, dr, rho, rhop, rhopp, uxc,   &
     &                      exc, ienergy, exmix)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ienergy
        integer, intent (in) :: ioption
        integer, intent (in) :: mesh

        real(kind=long), intent (in) :: dr
        real(kind=long), intent (in) :: exmix

        real(kind=long), intent (in), dimension (mesh) :: r
        real(kind=long), intent (in), dimension (mesh) :: rho
        real(kind=long), intent (in), dimension (mesh) :: rhop
        real(kind=long), intent (in), dimension (mesh) :: rhopp

! Output
        real(kind=long), intent (out) :: exc
        real(kind=long), intent (out), dimension (mesh) :: uxc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint

        real(kind=long) aln
        real(kind=long) ecp
        real(kind=long) ex
        real(kind=long) fx
        real(kind=long) fxc
        real(kind=long) rh
        real(kind=long) rs
        real(kind=long) rsx
        real(kind=long) zeta

        real(kind=long), dimension (2) :: cpot
        real(kind=long), dimension (2) :: d
        real(kind=long), dimension (:), allocatable :: dec
        real(kind=long), dimension (:), allocatable :: dex
        real(kind=long), dimension (:), allocatable :: dexc
        real(kind=long), dimension (2) :: dp
        real(kind=long), dimension (2) :: dpp
        real(kind=long), dimension (2) :: xpot

! Allocate Arrays
! ===========================================================================
        allocate (dec(mesh))
        allocate (dex(mesh))
        allocate (dexc(mesh))

! Procedure
! ===========================================================================
! Select which exchange-correlation potential is to be calculated
        select case (ioption)

! Wigner
!         case (1)
!          do ipoint = 1, mesh
!           rh = rho(ipoint)
!           call wigner (rh, ex, fx, exc, fxc)
!           uxc(ipoint) = fxc
!           dex(ipoint) = ex
!           dec(ipoint) = exc - ex
!          end do

! Hedin-Lundqvist
!         case (2)
!          dex = 0.d0
!          uxc = 0.d0
!          do ipoint = 1, mesh
!           rh = rho(ipoint)
!           if (rh .ne. 0.0d0) then
!            rs = 0.62035049d0*rh**(-0.3333333333333333d0)
!            rsx = rs/21.0d0
!            aln = dlog(1.0d0 + 1.0d0/rsx)
!            ecp = aln + (rsx**3*aln - rsx*rsx) + rsx/2 - 1.0d0/3.0d0
!            dex(ipoint) = -0.458175d0/rs - 0.0225d0*ecp
!            uxc(ipoint) = -0.6109d0/rs - 0.0225d0*aln
!           end if
!          end do

! Ceperley-Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
         case (3)
          do ipoint = 1, mesh
           rh = rho(ipoint)
           call cepal (rh, ex, fx, exc, fxc)
           uxc(ipoint) = fxc
           dex(ipoint) = ex
           dec(ipoint) = exc - ex
          end do

! Perdew exchange, Perdew correlation, generalized gradient approximation 1992
!         case (4)
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           dp(1) = 0.5d0*rhop(ipoint)
!           dpp(1) = 0.5d0*rhopp(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           dp(2) = 0.5d0*rhop(ipoint)
!           dpp(2) = 0.5d0*rhopp(ipoint)
!           call ggaxrad (3, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (2, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Becke exchange, Perdew correlation, generalized gradient approximation
!         case (5)
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           dp(1) = 0.5d0*rhop(ipoint)
!           dpp(1) = 0.5d0*rhopp(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           dp(2) = 0.5d0*rhop(ipoint)
!           dpp(2) = 0.5d0*rhopp(ipoint)
!           call ggaxrad (2, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (3, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Burke-Perdew-Ernzerhof generalized gradient approximation (1996)
!         case (6)
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           dp(1) = 0.5d0*rhop(ipoint)
!           dpp(1) = 0.5d0*rhopp(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           dp(2) = 0.5d0*rhop(ipoint)
!           dpp(2) = 0.5d0*rhopp(ipoint)
!           call ggaxrad (5, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (5, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Wigner-scaled LDA of PRA 46, R5320 (1992)
!         case (7)
!          do ipoint = 1, mesh
!           rh = rho(ipoint)
!           call wigscaled (rh, ex, fx, exc, fxc)
!           uxc(ipoint) = fxc
!           dex(ipoint) = ex
!           dec(ipoint) = exc - ex
!          end do

! Ceperley-Alder as parameterized by Perdew-Wang (1991)
!         case (8)
!          dp(1) = 0.0d0
!          dpp(1) = 0.0d0
!          dp(2) = 0.0d0
!          dpp(2) = 0.0d0
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           call ggaxrad (1, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (1, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Lee-Yang-Parr correlation
! Becke exchange, generalized gradient approximation
         case (9)
          do ipoint = 1, mesh
           d(1) = 0.5d0*rho(ipoint)
           dp(1) = 0.5d0*rhop(ipoint)
           dpp(1) = 0.5d0*rhopp(ipoint)
           d(2) = 0.5d0*rho(ipoint)
           dp(2) = 0.5d0*rhop(ipoint)
           dpp(2) = 0.5d0*rhopp(ipoint)
           call ggaxrad (2, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
           call ggacrad (4, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
           uxc(ipoint) = xpot(1) + cpot(1)
          end do

! Perdew-Wang exchange, generalized gradient approximation
!         case (10)
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           dp(1) = 0.5d0*rhop(ipoint)
!           dpp(1) = 0.5d0*rhopp(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           dp(2) = 0.5d0*rhop(ipoint)
!           dpp(2) = 0.5d0*rhopp(ipoint)
!           call ggaxrad (3, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (4, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Volko-Wilk-Nusair exchange-correlation (spin-polarization)
!         case (11)
!          do ipoint = 1, mesh
!           zeta = 0.0d0
!           d(1) = 0.5d0*rho(ipoint)*(1.0d0 + zeta)
!           d(2) = 0.5d0*rho(ipoint)*(1.0d0 - zeta)
!           call lsdavwn (d, dex(ipoint), dec(ipoint), xpot, cpot)
!           uxc(ipoint) = xpot(1) + cpot(1)
!          end do

! Lee-Yang-Parr correlation
! Becke exchange, generalized gradient approximation mixed with exact exchange.
!         case (12)
!          do ipoint = 1, mesh
!           d(1) = 0.5d0*rho(ipoint)
!           dp(1) = 0.5d0*rhop(ipoint)
!           dpp(1) = 0.5d0*rhopp(ipoint)
!           d(2) = 0.5d0*rho(ipoint)
!           dp(2) = 0.5d0*rhop(ipoint)
!           dpp(2) = 0.5d0*rhopp(ipoint)
!           call ggaxrad (2, r(ipoint), d, dp, dpp, xpot, dex(ipoint))
!           call ggacrad (4, r(ipoint), d, dp, dpp, cpot, dec(ipoint))
!           uxc(ipoint) = (1.0d0 - exmix)*xpot(1) + cpot(1)
!          end do
        end select

! ****************************************************************************
! Take care of the endpoints:
        uxc(1) = 2*uxc(2) - uxc(3)
        uxc(mesh) = 2*uxc(mesh-1) - uxc(mesh - 2)

! Calculate the total energy components
        if (ienergy .eq. 1) then
         if (ioption .ne. 12) then
          dexc = dex + dec
         else
          dexc = (1.0d0 - exmix)*dex + dec
         end if
         exc = 0.0d0
         do ipoint = 1, mesh
          exc = exc + dr*r(ipoint)**2*rho(ipoint)*(dexc(ipoint) - uxc(ipoint))
         end do
        end if

! Change the energies and potentials into Rydbergs - multiply by 2!
        exc = 2.0d0*exc
        uxc = 2.0d0*uxc

! Deallocate Arrays
! ===========================================================================
        deallocate (dec)
        deallocate (dex)
        deallocate (dexc)

! Format Statements
! ===========================================================================

        return
        end
