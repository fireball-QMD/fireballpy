! copyright info:
!
!                             @Copyright 1999
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
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth i
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
! factorial.f
! Program Description
! ===========================================================================
!
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function factorial (ifac)
        use precision
        implicit none
        real(kind=long) factorial
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: ifac
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer icount
        integer iprod
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (ifac .lt. 0) then
         stop ' Error in factorial  ----  ifac < 0 '
        else if (ifac .eq. 0 .or. ifac .eq. 1) then
         iprod = 1
        else
         iprod = 1
         do icount = 2, ifac
          iprod = iprod*icount
         end do
        end if
        factorial = real(iprod)
 
! Format Statements
! ===========================================================================
 
        return
        end