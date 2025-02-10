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

! clebsch_gordon.f90
! Program Description
! ===========================================================================
!       This routine calculates the Clebsch-Gordon coefficients which
! are represented by - <l1,l2;m1,m2|l,m>.
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
        function clebsch_gordon (l1, m1, l2, m2, l, m)
        use precision
        implicit none
        real(kind=long) clebsch_gordon

! Argument Declaration and Description
! ===========================================================================
        integer l, l1, l2
        integer m, m1, m2
!
! Local Parameters and Data Declaration
! ===========================================================================
! The maximum z value is to be only 6, since the lmax value is only 6 in the
! routine which calls this function.
        integer izmax
        parameter (izmax = 50)

! Local Variable Declaration and Description
! ===========================================================================
        integer iz

        real(kind=long) factorial
        real(kind=long) piece1
        real(kind=long) piece2
        real(kind=long) piece3
        real(kind=long) part3

        external factorial

! Procedure
! ===========================================================================
! Initialize the coefficient to zero.
        clebsch_gordon = 0.0d0

! First determine that the values of the l1, l2, and l satisfiy the
! triangle equation - |l1-l2| =< l >= l1 + l2.
        if (l .lt. abs(l1 - l2)) return
        if (l .gt. (l1 + l2)) return

! The other condition is that m = m1 + m2
        if (m .ne. (m1 + m2)) return

! The clebsch_gordon coefficient will be written as a product of three
! pieces. So clebsch_gordon = piece1*piece2*piece3
        piece1 = sqrt(((2.0d0*l + 1)*factorial(l1 + l2 - l)                 &
     &                *factorial(l1 - l2 + l)*factorial(-l1 + l2 + l))      &
     &                /factorial(l1 + l2 + l + 1))

        piece2 = sqrt(factorial(l1 + m1)*factorial(l1 - m1)                 &
     &                *factorial(l2 + m2)*factorial(l2 - m2)                &
     &                *factorial(l + m)*factorial(l - m))

        piece3 = 0.0d0
        do iz = 0, izmax
         if (((l1+l2-l-iz) .ge. -0.1) .and. ((l1-m1-iz) .ge. -0.1)          &
     &       .and. ((l2 + m2 - iz) .ge. -0.1)                               &
     &       .and. ((l - l2 + m1 + iz) .ge. -0.1)                           &
     &       .and. ((l - l1 - m2 + iz) .ge. -0.1)) then
          part3 = factorial(iz)*factorial(l1 + l2 - l - iz)                 &
     &             *factorial(l1 - m1 - iz)*factorial(l2 + m2 - iz)         &
     &             *factorial(l - l2 + m1 + iz)*factorial(l - l1 - m2 + iz)
          piece3 = piece3 + (-1)**iz/part3
         end if
        end do

! Now calculate the coefficient
        clebsch_gordon = piece1*piece2*piece3

! Format Statements
! ===========================================================================

        return
        end


        function factorial (ix)
        use precision
        real(kind=long) factorial

        integer ix
        integer k
        real(kind=long) fprod

        if (ix .lt. 0) then
         write(6,*) 'Stop !!!!   Factorial for negative x'
         stop
        end if
        fprod = 1.0
        do k = 1, ix
         fprod = fprod*real(k)
        end do
        factorial = fprod
        return
        end
