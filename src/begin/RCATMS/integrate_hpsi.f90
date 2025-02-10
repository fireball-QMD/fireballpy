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

! integrate_hpsi.f90
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
        subroutine integrate_hpsi (mesh, l, rcutoff, ebind, v, node, alnd)
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
        integer, intent (out) :: node

        real(kind=long), intent (out) :: alnd

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch
        integer ipoint
        integer istart

        real(kind=long) dr
        real(kind=long) h12
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
! Why is it multiplied by such an odd number such as 0.53 - is this supposed
! to be abohr?
        istart = 1
        imatch = mesh*0.53d0

! Initialize node counter
        node = 0

! Integrate out from origin
        u1 = 0.0d0
        u2 = dr
        v1 = 0.0d0
        v2 = v(istart+1) + v_ang(istart+1) + ebind
        do ipoint = istart + 2, imatch
         v3 = v(ipoint) + v_ang(ipoint) + ebind
         u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))             &
     &       /(1.0d0 - h12*v3)
         if (sign(1.0d0,u3) .ne. sign(1.0d0,u2)) node = node + 1
         u1 = u2
         u2 = u3
         v1 = v2
         v2 = v3
        end do
        v3 = v(imatch+1) + v_ang(imatch+1) + ebind
        u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)
        alnd = (u3 - u1)/(2.0d0*dr*u2)

! Now integrate in from rmax
        u1 = 0.0d0
        u2 = dr
        v1 = v(mesh) + v_ang(mesh) + ebind
        v2 = v(mesh-1) + v_ang(mesh-1) + ebind
        do ipoint = mesh-2, imatch, -1
         v3 = v(ipoint) + v_ang(ipoint) + ebind
         u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))             &
     &       /(1.0d0 - h12*v3)
         if (sign(1.0d0,u3) .ne. sign(1.0d0,u2)) node = node + 1
         u1 = u2
         u2 = u3
         v1 = v2
         v2 = v3
        end do
        v3 = v(mesh-imatch-1) + v_ang(mesh-imatch-1) + ebind
        u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)

! Calculate discontinuity in log derivative - it may look like this is lower
! order, but matching these actually matches the solution exactly
        alnd = alnd + (u3 - u1)/(2.0d0*dr*u2)

! Deallocate Arrays
! ===========================================================================
        deallocate (r)
        deallocate (v_ang)

! Format Statements
! ===========================================================================

        return
        end
