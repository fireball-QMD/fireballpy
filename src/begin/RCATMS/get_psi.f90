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

! get_psi.f90
! Program Description
! ===========================================================================
!       This routines integrates the Schroedinger equation for a trial binding
! energy. The number of nodes are returned and the discontinuity in the
! logarithmic derivative.
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
        subroutine get_psi (mesh, l, rcutoff, ebind, v, psi)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: mesh

        real(kind=long), intent (in) :: ebind
        real(kind=long), intent (in) :: rcutoff

        real(kind=long), intent (in), dimension (mesh) :: v

! Output
        real(kind=long), intent (out), dimension (mesh) :: psi

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch
        integer ipoint
        integer istart

        real(kind=long) anorm
        real(kind=long) dr
        real(kind=long) h12
        real(kind=long) ratio
        real(kind=long) u1, u2, u3
        real(kind=long) v1, v2, v3

        real(kind=long), dimension (:), allocatable :: r
        real(kind=long), dimension (:), allocatable :: v_ang

! Allocate Arrays
! ===========================================================================
        allocate (r(mesh))
        allocate (v_ang(mesh))

! Procedure
! ===========================================================================
        dr = rcutoff/(mesh - 1)
        h12 = dr**2/12.0d0
        do ipoint = 2, mesh
         r(ipoint) = (ipoint - 1)*dr
        end do
        v_ang(1) = 0.0d0
        v_ang(2:mesh) = l*(l + 1.0d0)/r(2:mesh)**2

! imatch = matching point for the integration
        istart = 1
        imatch = mesh/2

! Set up initial solution at origin and rmax
        psi(1) = 0.0d0
        psi(2) = dr
        psi(mesh) = 0.0d0
        psi(mesh-1) = dr

! Integrate out from origin
        u1 = 0.0d0
        u2 = dr
        v1 = v(istart) + ebind
        v2 = v(istart+1) + v_ang(istart+1) + ebind
        do ipoint = istart + 2, imatch
         v3 = v(ipoint) + v_ang(ipoint) + ebind
         u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.d0 - h12*v1))/(1.0d0 - h12*v3)
         psi(ipoint) = u3
         u1 = u2
         u2 = u3
         v1 = v2
         v2 = v3
        end do
        ratio = u3

! Now integrate in from rmax
        u1 = 0.0d0
        u2 = dr
        v1 = v(mesh) + v_ang(mesh) + ebind
        v2 = v(mesh-1) + v_ang(mesh-1) + ebind
        do ipoint = mesh-2, imatch, -1
         v3 = v(ipoint) + v_ang(ipoint) + ebind
         u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)
         psi(ipoint) = u3
         u1 = u2
         u2 = u3
         v1 = v2
         v2 = v3
        end do
        ratio = ratio/u3

! Match solutions
        psi(imatch:mesh) = psi(imatch:mesh)*ratio

! normalize using trapezoidal rule -- remember the end points are zero
        anorm = 1.0d0/sqrt(dot_product(psi,psi)*dr)
        psi = anorm*psi

! Deallocate Arrays
! ===========================================================================
        deallocate (r)
        deallocate (v_ang)

! Format Statements
! ===========================================================================

        return
        end
