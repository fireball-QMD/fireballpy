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
! ggaxrad.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
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
        subroutine ggaxrad (mode, rin, rho, rhop, rhopp, xpot, xen)
        use constants
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode
        real(kind=long), intent (in) :: rin
        real(kind=long), intent (in), dimension (2) :: rho
        real(kind=long), intent (in), dimension (2) :: rhop
        real(kind=long), intent (in), dimension (2) :: rhopp
! Output
        real(kind=long), intent (out) :: xen
        real(kind=long), intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin
        real(kind=long) density
        real(kind=long) densityp
        real(kind=long) densitypp
        real(kind=long) ex
        real(kind=long) fermik
        real(kind=long) r
        real(kind=long) s
        real(kind=long) u
        real(kind=long) v
        real(kind=long) vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
         if (rho(ispin) .le. eps) then
          xpot(ispin) = 0.0d0
         else
          density = 2.0d0*rho(ispin)
          if (mode .eq. 1) then
           call xlda (density, vx, ex)
          else if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 5) then
           densityp = 2.0d0*rhop(ispin)
           densitypp = 2.0d0*rhopp(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)
           s = abs(densityp)/(fermik*density)
           u = abs(densityp)*densitypp/(density*density*fermik**3)
           v = (densitypp + 2.0d0*densityp/r)/(density*fermik*fermik)
           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho(ispin)*ex
         end if
        end do

! Energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================

        return
        end
