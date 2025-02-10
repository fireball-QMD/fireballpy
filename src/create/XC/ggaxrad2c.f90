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
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!
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


! ggaxrad2c.f90
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
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggaxrad2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                        rhopz, xpot, xen)
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
        real(kind=long), intent (in), dimension (2) :: rhopz
        real(kind=long), intent (in), dimension (2) :: rhoz
        real(kind=long), intent (in), dimension (2) :: rhozz

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
        real(kind=long) densitypz
        real(kind=long) densityz
        real(kind=long) densityzz
        real(kind=long) ex
        real(kind=long) fermik
        real(kind=long) pi
        real(kind=long) r
        real(kind=long) s
        real(kind=long) u
        real(kind=long) v
        real(kind=long) vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        pi = 3.141592653589793238462643D0

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
           densityz = 2.0d0*rhoz(ispin)
           densityzz = 2.0d0*rhozz(ispin)
           densitypz = 2.0d0*rhopz(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
           s = sqrt(densityp**2 + densityz**2)/(fermik*density)
           u = (densityp**2*densitypp + 2.0d0*densityp*densityz*densitypz    &
     &          + densityz**2*densityzz)/(s*density**3*fermik**4)
           v = (densitypp + densityp/r + densityzz)/(density*fermik*fermik)

           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad2c : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho(ispin)*ex
         end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end