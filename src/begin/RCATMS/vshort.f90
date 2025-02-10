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

! vshort.f90
! Program Description
! ===========================================================================
!       This function returns the value of the core potential in real space.
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
        function vshort (r)
        use precision
        use pp_storage
        implicit none
        real(kind=long) vshort

! Argument Declaration and Description
! ===========================================================================
        real(kind=long) r

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: norder = 3

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid, iright
        integer iprod
        integer isum

        real(kind=long) prod

! Procedure
! ===========================================================================
! Initialize vshort to zero
        vshort = 0.0d0

! Now the shortranged local potential
! If r is slightly negative return the value at r = 0
        if (r .lt. 0.0d0) then
         vshort = v_short(1)
         return
        end if

! The short stuff is zero at long range.
        if (r .gt. rrc_short) return

        imid = int(r/drr_short) + 1

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
        else if (iright .gt. npoints_short) then
         ileft = npoints_short - (norder+1)
         iright = npoints_short
        end if

! Now interpolate
        do isum = ileft, iright
         prod = 1.0d0
         do iprod = ileft, iright
          if (iprod .ne. isum) then
           prod = prod*(r - r_short(iprod))/(r_short(isum) - r_short(iprod))
          end if
         end do
         vshort = vshort + v_short(isum)*prod
        end do

! Format Statements
! ===========================================================================

        return
        end
