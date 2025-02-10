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

! psirc.f90
! Program Description
! ===========================================================================
!       This routine finds the wavefunction and energy for a potential given
! the boundary condition that psi = 0 at rcutoff. (Hartree, Bohr units). This
! routine can do any l state. The eigenvalue is negative for a bound state.
! This routine will solve for any excited state, use nexcit = 0 for the usual
! node_loweress case. Non-zero nexcit will give the corresponding excited state
! with angular momentum l.
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
        subroutine psirc (mesh, nexcite, l, rcutoff, v, eout, psi)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: mesh
        integer, intent (in) :: nexcite

        real(kind=long), intent (in) :: rcutoff

        real(kind=long), intent (in), dimension (mesh) :: v

! Output
        real(kind=long), intent (out) :: eout

        real(kind=long), intent (out), dimension (mesh) :: psi

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: small = 1.0d-14

! Local Variable Declaration and Description
! ===========================================================================
        integer iteration
        integer node
        integer node_lower, node_upper

        real(kind=long) alnd, alnd1, alnd2
        real(kind=long) alnd_lower, alnd_upper
        real(kind=long) ebind, ebind1, ebind2
        real(kind=long) ebind_lower, ebind_upper

        logical test

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Assume all energies are not orders of magnitude away from atomic units --
! otherwise you might want to put a guessed energy here to start.
! Find upper bound to binding energy
        ebind_upper = 1.0d0
        do iteration = 1, 100
         ebind = 2.0d0*ebind_upper
         call integrate_hpsi (mesh, l, rcutoff, ebind, v, node, alnd)
         if (test(node,nexcite,alnd)) exit
         ebind_upper = ebind
        end do
        if (.not. test(node,nexcite,alnd))                                  &
     &   stop ' no upper bound found -- stopping '
        node_upper = node
        alnd_upper = alnd
        ebind_upper = ebind

! Find lower bound to binding energy
        ebind_lower = -1.0d0
        do iteration = 1, 100
         ebind = 2.0d0*ebind_lower
         call integrate_hpsi (mesh, l, rcutoff, ebind, v, node, alnd)
         if (.not. test(node,nexcite,alnd)) exit
         ebind_lower = ebind
        end do
        if (test(node,nexcite,alnd)) stop ' no upper bound found -- stopping '
        node_lower = node
        alnd_lower = alnd
        ebind_lower = ebind

! Do binary chop to get close
        do iteration = 1, 100
         if (node_lower .eq. node_upper) exit
         ebind = 0.5d0*(ebind_upper + ebind_lower)
         call integrate_hpsi (mesh, l, rcutoff, ebind, v, node, alnd)
         if (test(node,nexcite,alnd)) then
          ebind_upper = ebind
          node_upper = node
          alnd_upper = alnd
         else
          ebind_lower = ebind
          node_lower = node
          alnd_lower = alnd
         end if
        end do
        if (node_lower .ne. node_upper) then
         write (*,*) ' ******************* WARNING ************************* '
         write (*,*) ' No convergence in binary chop. Generally, this means  '
         write (*,*) ' your grid is too course or something else is causing  '
         write (*,*) ' an extra node. '
        end if

! Linearly interpolate to get good eigenvalues
        ebind1 = ebind_upper
        ebind2 = ebind_lower
        alnd1 = alnd_upper
        alnd2 = alnd_lower
        do iteration = 1, 100
         if (abs(alnd1-alnd2) .lt. small) exit
         ebind = ebind1 - (ebind2 - ebind1)*alnd1/(alnd2 - alnd1)
         if (ebind .gt. ebind_upper .or. ebind .lt. ebind_lower)          &
     &    ebind = 0.5d0*(ebind_upper + ebind_lower)
         call integrate_hpsi (mesh, l, rcutoff, ebind, v, node, alnd)
         if (test(node,nexcite,alnd)) then
          ebind_upper = ebind
          alnd_upper = alnd
         else
          ebind_lower = ebind
          alnd_lower = alnd
         end if
         alnd2 = alnd1
         alnd1 = alnd
         ebind2 = ebind1
         ebind1 = ebind
        end do
        if (abs(alnd1-alnd2) .gt. small)                                  &
     &   stop ' No convergence in linear extrapolation '
        eout = -ebind
        call get_psi (mesh, l, rcutoff, ebind, v, psi)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end


        function test (node, nexcite, alnd)
        use precision
        implicit none
        logical test

        integer, intent (in) :: node
        integer, intent (in) :: nexcite

        real(kind=long), intent (in) :: alnd

        test = (node .lt. nexcite) .or.                                    &
     &         (node .eq. nexcite .and. alnd .gt. 0.0d0)

        return
        end
