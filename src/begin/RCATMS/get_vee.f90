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

! get_vee.f90
! Program Description
! ===========================================================================
!
!       This subroutine evaluates the electron-electron Hartree potential for
! a given density. It is given by:
!
!    vee(r) = 2*int(dr' n(r)/!r-r'!)
!           = 2*(1/r*int[0,r] sigma(t)dt + int[r,rc] sigma(t)/t dt  )
!
! Let x1(r) = int[0,r] sigma(t)dt
!     x2(r) = int[r,rc] sigma(t)/t dt.
!
! then, vee(r) = 2/r*x1(r) + 2*(x2(mesh) - x2(i))
!
! Integrate to get x1(r), via trapezoidal rule.
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
        subroutine get_vee (mesh, dr, sigma, vee)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mesh

        real(kind=long), intent (in) :: dr
        real(kind=long), intent (in), dimension (mesh) :: sigma

! Output
        real(kind=long), intent (out), dimension (mesh) :: vee

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint

        real(kind=long) a1, a2
        real(kind=long) b1, b2
        real(kind=long) r

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Find x1(r)
        a1 = 0.0d0
        vee = 0.0d0
        do ipoint = 2, mesh
         r = (ipoint - 1)*dr
         a2 = 0.5d0*sigma(ipoint)
         a1 = a1 + a2
         vee(ipoint) = 2.0d0*a1*dr/r
         a1 = a1 + a2
        end do

! Next get x2(r) - also trapezoidal rule
        b1 = 0.0d0
        do ipoint = mesh, 2, -1
         r = (ipoint - 1)*dr
         b2 = 0.5d0*sigma(ipoint)/r
         b1 = b1 + b2
         vee(ipoint) = vee(ipoint) + 2.0d0*b1*dr
         b1 = b1 + b2
        end do
        vee(1) = vee(1) + 2.0d0*b1*dr

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
