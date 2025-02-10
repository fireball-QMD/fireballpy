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

! psiofr.f90
! Program Description
! ===========================================================================
!       This function returns the values psiofr(r) for the corresponding
! shell.  The radial functions are normalized as:

!  int ( psiofr**2  r**2  dr ) = 1.0

! The wavefunctions for each atom must be read in by calling readpsi
! separately for each atom type.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Code rewritten by:
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
        function psiofr (issh, r, nssh, mesh, npoints, drr, rrc, rr, psi)
        use precision
        implicit none
        real(kind=long) psiofr

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: issh
        integer, intent(in) :: mesh
        integer, intent(in) :: nssh

        integer, intent(in), dimension (nssh) :: npoints

        real(kind=long), intent (in) :: r

        real(kind=long), intent(in), dimension (nssh) :: drr
        real(kind=long), intent(in), dimension (nssh, mesh) :: psi
        real(kind=long), intent(in), dimension (nssh, mesh) :: rr
        real(kind=long), intent(in), dimension (nssh) :: rrc

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: norder = 3

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid, iright
        integer iprod
        integer isum

        real(kind=long) prod
        real(kind=long) sum

! Procedure
! ===========================================================================
! Initialize to zero
        psiofr = 0.0d0
        if (r .gt. rrc(issh)*0.529177249d0) return

! If r is slightly negative return the value at r = 0
        if (r .lt. 0.0d0) then
         psiofr = psi(issh,1)
         return
        end if

        imid = idint(r/drr(issh)) + 1

! Find starting and ending points for the interpolation
        if (mod(norder + 1, 2) .eq. 0) then
         ileft = imid - ((norder-1)/2)
         iright = imid + ((norder+1)/2)
        else
         ileft = imid - (norder/2)
         iright = imid + (norder/2)
        end if

        if (ileft .lt. 1) then
         ileft = 1
         iright = norder + 1
        else if (iright .gt. npoints(issh)) then
         ileft = npoints(issh) - (norder+1)
         iright = npoints(issh)
        end if

! Now interpolate
        sum = 0.0d0
        do isum = ileft, iright
         prod = 1.0d0
         do iprod = ileft, iright
          if (iprod .ne. isum) then
           prod = prod*(r - rr(issh,iprod))/(rr(issh,isum) - rr(issh,iprod))
          end if
         end do
         sum = sum + psi(issh,isum)*prod
        end do

        psiofr = sum

! Format Statements
! ===========================================================================

        return
        end
