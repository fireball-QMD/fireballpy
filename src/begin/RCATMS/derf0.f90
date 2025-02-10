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

! derf0.f90
! Program Description
! ===========================================================================
!       Error function, see Abramowitz and Stegun.
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
        function derf0 (x)
        use precision
        implicit none
        real(kind=long) derf0

! Argument Declaration and Description
! ===========================================================================
        real(kind=long) x

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real(kind=long) a1, a2, a3, a4, a5
        real(kind=long) p
        real(kind=long) t
        real(kind=long) y, z

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
        y = x
        if (y .lt. 0.0d0) y = -y
        if (y .gt. 1.0d-7) then
         p = 0.3275911d0
         t = 1.0d0/(1.0d0 + p*y)
         a1 = 0.254829592d0
         a2 = -0.284496736d0
         a3 = 1.421413741d0
         a4 = -1.453152027d0
         a5 = 1.061405429d0
         z = 1.0d0 - exp(-y*y)*(a1*t + a2*t*t + a3*(t**3) + a4*(t**4)       &
     &                                                    + a5*(t**5))
         if (x .lt. 0.0d0) z = -z
        else
         z = 0.0
        end if
        derf0 = z


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
