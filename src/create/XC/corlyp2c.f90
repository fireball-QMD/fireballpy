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


! corlyp2c.f90
! Program Description
! ========================================================================
!
!  The subroutine corlyp calculates the correlation energy using
!  the Lee, Yang and Parr generalized gradient approximation.
!
!
! Input Variables
!
! tpot .... T  evaluate correlation energy and potential
!           F  evaluate energy only (a posteriori scheme)
! x ....... dummy
! pa ...... spin up density
! pb ...... spin down density
! dpaofr .. 1st partial derivative of spin up with respect to r
! dpaofz .. 1st partial derivative of spin up with respect to z
! d2paofr . 2nd partial derivative of spin up with respect to r
! d2paofz . 2nd partial derivative of spin up with respect to z
! dpbofr .. 1st partial derivative of spin down with respect to r
! dpbofz .. 1st partial derivative of spin down with respect to z
! d2pbofr . 2nd partial derivative of spin down with respect to r
! d2pbofz . 2nd partial derivative of spin down with respect to z
!
! Output Variables
!
! ec ...... correlation energy per electron
! vp ...... correlation potential
!
! ===========================================================================
! Code written by:
! Richard B. Evans and James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine corlyp2c (tpot, r, pa, pb, dpaofr, dpbofr, d2paofr,     &
     &                       d2pbofr, dpaofz, dpbofz, d2paofz, d2pbofz,    &
     &                       ec, vp)
        use precision, only: wp
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        logical  tpot

        real(kind=wp) r
        real(kind=wp) pa
        real(kind=wp) pb
        real(kind=wp) dpaofr
        real(kind=wp) dpaofz
        real(kind=wp) dpbofr
        real(kind=wp) dpbofz
        real(kind=wp) d2paofr
        real(kind=wp) d2paofz
        real(kind=wp) d2pbofr
        real(kind=wp) d2pbofz

! Output
        real(kind=wp) ec
        real(kind=wp) vp

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=wp), parameter :: aa = 0.04918d0
        real(kind=wp), parameter :: bb = 0.132d0
        real(kind=wp), parameter :: cc = 0.2533d0
        real(kind=wp), parameter :: dd = 0.349d0
        real(kind=wp), parameter :: t13 = 1.0d0/3.0d0
        real(kind=wp), parameter :: t23 = 2.0d0/3.0d0
        real(kind=wp), parameter :: t43 = 4.0d0/3.0d0
        real(kind=wp), parameter :: t53 = 5.0d0/3.0d0
        real(kind=wp), parameter :: t73 = 7.0d0/3.0d0
        real(kind=wp), parameter :: t83 = 8.0d0/3.0d0

! Local Variable Delcaration and Descrition
! ======================================================================
! constants
        real(kind=wp) Cf
        real(kind=wp) pi
        real(kind=wp) dep
        real(kind=wp) dp_one
        real(kind=wp) ep

! density and derivatives
        real(kind=wp) p             ! Sum of spin up and spin down densities
        real(kind=wp) dpr           ! Sum of spin up and down partials wrt r
        real(kind=wp) dpz           ! Sum of spin up and down partials wrt z
        real(kind=wp) d2pr          ! Sum of spin 2nd partials wrt r
        real(kind=wp) d2pz          ! Sum of spin 2nd partials wrt z

        real(kind=wp) gamma
        real(kind=wp) dgammap       ! partial derivative of gamma wrt pa
        real(kind=wp) dgammar       ! partial derivative of gamma wrt r
        real(kind=wp) dgammaz       ! partial derivative of gamma wrt z
        real(kind=wp) d2gammaz      ! 2nd Partial derivative of gamma wrt z
        real(kind=wp) d2gammar      ! 2nd Partial derivative of gamma wrt r

        real(kind=wp) Fofp
        real(kind=wp) dFofp         ! partial derivative of Fofp wrt pa
        real(kind=wp) dFofpz        ! partial derivative of Fofp wrt z
        real(kind=wp) dFofpr        ! partial derivative of Fofp wrt r
        real(kind=wp) d2Fofpr       ! 2nd Partial derivative of F wrt r
        real(kind=wp) d2Fofpz       ! 2nd Partial derivative of F wrt z

        real(kind=wp) Gofp
        real(kind=wp) dGofp         ! partial derivative of Gofp wrt pa
        real(kind=wp) dGofpz        ! partial of Gofp wrt z
        real(kind=wp) dGofpr        ! partial of Gofp wrt r
        real(kind=wp) d2Gofpr       ! 2nd partial of Gofp wrt r
        real(kind=wp) d2Gofpz       ! 2nd partial of Gofp wrt z
        real(kind=wp) d2Gofp        ! Laplacian of Gofp

        real(kind=wp) tW            ! Weizsacker kinetic-energy density for p
        real(kind=wp) taW           ! Weizsacker K.E. dens. for spin up
        real(kind=wp) tbW           ! Weizsacker K.E. dens. for spin down

! Procedure
! ===========================================================================
! The notation used for our variables is similar to the original
! notation used by Lee, Yang and Parr. (Physical Review B 37 785 (1988))
! Initialize pi
        pi = 3.141592653589793238462643D0

        Cf = 0.3d0*(3.0d0*pi**2)**t23

! Here the option to calculate the potential is evaluated
        if (tpot) then
         p = pa + pb

         dpr = dpaofr + dpbofr
         dpz = dpaofz + dpbofz

         d2pr = d2paofr + d2pbofr
         d2pz = d2paofz + d2pbofz

         ep = p**(-t53)*exp(-cc*p**(-t13))
         dep = cc*t13*p**(-t43) - t53/p

         dp_one = 1.0d0 + dd*p**(-t13)

         gamma = 2.0d0*(1.0d0 - (pa**2 + pb**2)/p**2)
         dgammap = - 4.0d0*pa/p**2 - 2.0d0*(gamma - 2.0d0)/p
         dgammar = - 4.0d0*(pa*dpaofr + pb*dpbofr)/p**2                     &
     &             - 2.0d0*dpr*(gamma - 2.0d0)/p
         dgammaz = - 4.0d0*(pa*dpaofz + pb*dpbofz)/p**2                     &
     &             - 2.0d0*dpz*(gamma - 2.0d0)/p
         d2gammar =                                                         &
     &    - 4.0d0*(pa*d2paofr + dpaofr**2 + pb*d2pbofr + dpbofr**2)/p**2    &
     &    + 8.0d0*(pa*dpaofr + pb*dpbofr)*dpr/p**3                          &
     &    - 2.0d0*d2pr*(gamma - 2.0d0)/p                                    &
     &    - 2.0d0*dpr*(p*dgammar - (gamma - 2.0d0)*dpr)/p**2
         d2gammaz =                                                         &
     &    - 4.0d0*(pa*d2paofz + dpaofz**2 + pb*d2pbofz + dpbofz**2)/p**2    &
     &    + 8.0d0*(pa*dpaofz + pb*dpbofz)*dpz/p**3                          &
     &    - 2.0d0*d2pz*(gamma - 2.0d0)/p                                    &
     &    - 2.0d0*dpz*(p*dgammaz - (gamma - 2.0d0)*dpz)/p**2

         Fofp  = gamma/dp_one
         dFofp = (dgammap + dd*t13*p**(-t43)*Fofp)/dp_one
         dFofpr = (dgammar + dd*t13*p**(-t43)*Fofp*dpr)/dp_one
         dFofpz = (dgammaz + dd*t13*p**(-t43)*Fofp*dpz)/dp_one
         d2Fofpr =                                                          &
     &    (dp_one*(d2gammar - dd*t13*t43*p**(-t73)*Fofp*dpr**2              &
     &             + dd*t13*p**(-t43)*(dFofpr*dpr + Fofp*d2pr))             &
     &     + dd*t13*p**(-t43)*dpr*(dgammar + dd*t13*p**(-t43)*Fofp*dpr))    &
     &    /dp_one**2
         d2Fofpz =                                                          &
     &    (dp_one*(d2gammaz - dd*t13*t43*p**(-t73)*Fofp*dpz**2              &
     &             + dd*t13*p**(-t43)*(dFofpz*dpz + Fofp*d2pz))             &
     &    + dd*t13*p**(-t43)*dpz*(dgammaz + dd*t13*p**(-t43)*Fofp*dpz))     &
     &    /dp_one**2

         Gofp  = Fofp*ep
         dGofp = dFofp*ep + Gofp*dep
         dGofpr = dFofpr*ep + Gofp*dep*dpr
         dGofpz = dFofpz*ep + Gofp*dep*dpz
         d2Gofpr =                                                          &
     &    (d2Fofpr + dFofpr*dpr*dep)*ep + (dGofpr*dpr + Gofp*d2pr)*dep      &
     &     + Gofp*dpr**2*(t53/p**2 - cc*t13*t43*p**(-t73))
         d2Gofpz =                                                          &
     &    (d2Fofpz + dFofpz*dpz*dep)*ep + (dGofpz*dpz + Gofp*d2pz)*dep      &
     &     + Gofp*dpz**2*(t53/p**2 - cc*t13*t43*p**(-t73))

         d2Gofp = d2Gofpr + dGofpr/r + d2Gofpz

         vp = - aa*(dFofp*p + Fofp)                                         &
     &    - 2.0d0**t53*aa*bb*Cf                                             &
     &                *(dGofp*(pa**t83 + pb**t83) + t83*Gofp*pa**t53)       &
     &    - aa*bb/4.0d0*(p*d2Gofp + 4.0d0*(dGofpr*dpr + dGofpz*dpz)         &
     &                   + 4.0d0*Gofp*(d2pr + dpr/r + d2pz)                 &
     &                   + dGofp*(p*(d2pr + dpr/r + d2pz)                   &
     &                   - (dpr**2 + dpz**2)))                              &
     &    - aa*bb/36.0d0                                                    &
     &        *(3.0d0*pa*d2Gofp + 4.0d0*(dpaofr*dGofpr + dpaofz*dGofpz)     &
     &          + 4.0d0*Gofp*(d2paofr + dpaofr/r + d2paofz)                 &
     &          + 3.0d0*dGofp*(pa*(d2paofr + dpaofr/r + d2paofz)            &
     &                         + pb*(d2pbofr + dpbofr/r + d2pbofz))         &
     &          + dGofp*(dpaofr**2 + dpaofz**2 + dpbofr**2 + dpbofz**2))
        else
         vp = 0.0d0
        end if

! correlation energy per electron
        tW = ((dpr**2 + dpz**2)/p - (d2pr + dpr/r + d2pz))/8.0d0
        taW = ((dpaofr**2 + dpaofz**2)/pa                                   &
     &         - (d2paofr + dpaofr/r + d2paofz))/8.0d0
        tbW = ((dpbofr**2 + dpbofz**2)/pb                                   &
     &         - (d2pbofr + dpbofr/r + d2pbofz))/8.0d0

        ec = 2.0d0**t23*Cf*(pa**t83 + pb**t83) - p*tW                       &
     &      + (pa*taW + pb*tbW)/9.0d0                                       &
     &      + (pa*(d2paofr + dpaofr/r + d2paofz)                            &
     &         + pb*(d2pbofr + dpbofr/r + d2pbofz))/18.0d0
        ec = - aa*gamma/dp_one*(p + 2.0d0*bb*p**(-t53)*exp(-cc*p**(-t13))*ec)
        ec = ec/p

! Format Statements
! ===========================================================================

        return
        end
