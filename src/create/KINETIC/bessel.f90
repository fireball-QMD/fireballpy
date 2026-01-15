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
!                      Kirk VanOpdorp
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
!
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

 
! bessel.f
! Program Description
! ===========================================================================
!
!       spherical bessel function
!
! ===========================================================================
! Original code written by Otto F. Sankey
 
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real(kind=wp) function jl2 (l,x)
        use precision, only: wp
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer l
 
        real(kind=wp) x
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (l .gt. 6) then
         write (*,*) ' You are trying to calculate a spherical bessel '
         write (*,*) ' function for L = 7 or higher.  Currently, terms '
         write (*,*) ' for L > 6 are not considered in this function. '
         stop 'error in bessel'
        end if
 
! L = 0
! j0 = sinx/x
        if (l .eq. 0) then
         if (x .gt. 1.0d-4) then
          jl2 = dsin(x)/x
         else
          jl2 = (1.0d0 - x**2/6.0d0 + x**4/120.0d0)
         end if
        end if
 
! L = 1
! j1 = sinx/x**2 - cosx/x
        if (l .eq. 1) then
         if (x .gt. 1.0d-1) then
          jl2 = dsin(x)/(x**2) - cos(x)/x
         else
          jl2 = 1.0d0 - x**2/10.0d0 + x**4/280.0d0 - x**6/15120.0d0 &
     &                + x**8/1330560.0d0
          jl2 = (x/3.0d0)*jl2
         end if
        end if
 
! L = 2
! j2 = (3/x**3 - 1/x)*sinx - 3*cosx/x**2
        if (l .eq. 2) then
         if (x .gt. 5.0d-1) then
          jl2 = (3.0d0/x**3 - 1.0d0/x)*dsin(x) - (3.0d0/x**2)*dcos(x)
         else
          jl2 = 1.0d0 - x**2/14.0d0 + x**4/504.0d0 - x**6/33264.0d0 &
     &                + x**8/3459456.0d0
          jl2 = (x**2/15.0d0)*jl2
         end if
        end if
 
! L = 3
        if (l .eq. 3) then
         if (x .gt. 5.0d-1) then
           jl2 = (15.0d0/x**4 - 6.0d0/x**2)*dsin(x) &
     &         - (15.0d0/x**3 - 1.0d0/x)*dcos(x)
         else
           jl2 = 1.0d0 - x**2/18.0d0 + x**4/792.0d0 - x**6/61776.0d0 &
     &                + x**8/7413120.0d0
           jl2 = (x**3/105.0d0)*jl2
         end if
        end if
 
! L = 4
        if (l .eq. 4) then
         if (x .gt. 5.0d-1) then
           jl2 = (105.0d0/x**5 - 45.0d0/x**3 + 1.0d0/x)*dsin(x) &
     &         - (105.0d0/x**4 - 10.0d0/x**2)*dcos(x)
         else
           jl2 = 1.0d0 - x**2/22.0d0 + x**4/1144.0d0 - x**6/102960.0d0 &
     &                + x**8/14002560.0d0
           jl2 = (x**4/945.0d0)*jl2
         end if
        end if
 
! L = 5
        if (l .eq. 5) then
         if (x .gt. 5.0d-1) then
          jl2 = (945.0d0/x**6 - 420.0d0/x**4 + 15.0d0/x**2)*dsin(x) &
     &         - (945.0d0/x**5 - 105.0d0/x**3 + 1.0d0/x)*dcos(x)
         else
          jl2 = 1.0d0 - x**2/26.0d0 + x**4/1560.0d0 - x**6/159120.0d0
          jl2 = (x**5/10395.0d0)*jl2
         end if
        end if
 
! L = 6
        if (l .eq. 6) then
         if (x .gt. 5.0d-1) then
          jl2 = (10395.0d0/x**7 - 4725.0d0/x**5 + 210.0d0/x**3 &
     &                                          - 1.0d0/x)*dsin(x) &
     &         - (10395.0d0/x**6 - 1260.0d0/x**4 + 21.0d0/x**2)*dcos(x)
         else
          jl2 = 1.0d0 - x**2/30.0d0 + x**4/2040.0d0 - x**6/232560.0d0
          jl2 = (x**6/135135.0d0)*jl2
         end if
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end