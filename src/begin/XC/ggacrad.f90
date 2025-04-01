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

! ggacrad.f90
! Program Description
! ===========================================================================
!      This routine calculates the correlation potential and energy density.
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
        subroutine ggacrad (mode, rin, rho, rhop, rhopp, cpot, cen)
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
        real(kind=long), intent (out) :: cen

        real(kind=long), intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: crs = 1.91915829267751281d0
        real(kind=long), parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        real(kind=long) alfc
        real(kind=long) density
        real(kind=long) densityp
        real(kind=long) densityp11
        real(kind=long) densityp12
        real(kind=long) densityp22
        real(kind=long) densitypp
        real(kind=long) ec
        real(kind=long) ecrs
        real(kind=long) eczet
        real(kind=long) fermik
        real(kind=long) g
        real(kind=long) gsfermik
        real(kind=long) h
        real(kind=long) r
        real(kind=long) rs
        real(kind=long) sfermik
        real(kind=long) t
        real(kind=long) uu
        real(kind=long) vv
        real(kind=long) ww
        real(kind=long) zet
        real(kind=long) ztp

        real(kind=long), dimension (2) :: dvc, vc
        real(kind=long), dimension (2) :: flip

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! LSDA
        density = rho(1) + rho(2)
        densityp = rhop(1) + rhop(2)
        densitypp = rhopp(1) + rhopp(2)
        cen = 0.0d0
        if (density .le. eps) then
         cen = 0.0d0
         cpot(1) = 0.0d0
         cpot(2) = 0.0d0
         return
        end if

        if (mode .ne. 4) then
         zet = (rho(1) - rho(2))/density
         ztp = (rhop(1) - rhop(2) - zet*densityp)/density
         fermik = (3.0d0*pi*pi*density)**(1.0d0/3.0d0)
         rs = crs/fermik
         call corlsd (rs, zet, ec, vc(1), vc(2), ecrs, eczet, alfc)

! GGA correction to LSDA
         select case (mode)
          case (1)
           dvc = 0.0d0
           h = 0.0d0
          case (2)
            stop
   !        sfermik = 2.0d0*sqrt(fermik/pi)
   !        g = ((1.0d0 + zet)**(2.0d0/3.0d0)                               &
   !  &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
   !        gsfermik = 2.0d0*sfermik*g
   !        t = abs(densityp)/(density*gsfermik)
   !        uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
   !        vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
   !        ww = densityp*ztp/(density*gsfermik*gsfermik)
   !        call corgga (rs, zet, t, uu, vv, ww, h, dvc(1), dvc(2))
          case (3)
           uu = abs(densityp)*densitypp
           vv = densitypp + 2.0d0*densityp/r
           densityp11 = rhop(1)*rhop(1)
           densityp22 = rhop(2)*rhop(2)
           densityp12 = rhop(1)*rhop(2)
           if (rhop(1) .ne. 0.0d0 .or. rhop(2) .ne. 0.0d0)then
            call corga86 (rho(1),rho(2), densityp11, densityp22,           &
     &                    densityp12, uu, vv, h, dvc(1), dvc(2))
           end if
          case (5)
           sfermik = 2.0d0*sqrt(fermik/pi)
           g = ((1.0d0 + zet)**(2.0d0/3.0d0)                               &
     &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
           gsfermik = 2.0d0*sfermik*g
           t = abs(densityp)/(density*gsfermik)
           uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
           vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
           ww = densityp*ztp/(density*gsfermik*gsfermik)
           call corpbe (rs, zet, t, uu, vv, ww, 1, 1, ec, vc(1), vc(2), h,   &
     &                  dvc(1),dvc(2))
          end select
          cpot = vc + dvc
          cen = ec + h
        else if (mode .eq. 4) then
         flip = rhopp + 2.0d0*rhop/r
         call corlyp (.true., rho(1), rho(2), rhop(1), rhop(2), flip(1),     &
     &                flip(2), cen, cpot(1), cpot(2))
         !write (*,*) 'CONA', cen, cpot
        else
         stop 'ggacrad : mode improper'
        end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
