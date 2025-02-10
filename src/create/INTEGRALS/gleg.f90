! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

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


! gleg.f90
! Program Description
! ===========================================================================
!       Calculates the roots of the Legendre polynomials.
!
! ===========================================================================
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
        subroutine gleg (ctheta, ctheta_weights, ntheta_max)
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ntheta_max

! Output
        real(kind=long), intent (out), dimension (ntheta_max) :: ctheta 
        real(kind=long), intent (out), dimension (ntheta_max) :: ctheta_weights 

! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iterm
        integer itheta
        integer jterm
        integer nw
        integer ndown, nup
        integer node

        real(kind=long) pl1
        real(kind=long) pl2
        real(kind=long) pl3
        real(kind=long) plp
        real(kind=long) pdown, pup
        real(kind=long) xx

! Procedure
! ===========================================================================
! Calculate the positions and weights for an n point Gauss-Legendre
! integration formula. The positions are returned in x and the weights in w.
! The integration interval is from -1 to 1.
!
! Use symmetry -- only need 1/2 the points
! Find the nodes of P sub n.
        xx = 1.0d0
        nw = (ntheta_max + 1)/2
        pl2 = 2.0d0
        plp = 1.0d0

        ndown = nw
        nup = 0
        do iterm = 1, nw
 
! The variables pup and pdown are upper and lower bounds for the nodes.
         pup = xx
         pdown = 0.0d0
         do jterm = 1, 50
 
! Use either binary chop or Newton-Raphson to get to node. Use the property 
! of the recursion relation that the number of sign changes tells how many 
! nodes are located between xx and 1.
          if (nup .eq. iterm - 1 .and. ndown .eq. iterm) then
           xx = xx - pl2/plp
           if (xx .lt. pdown .or. xx .gt. pup) xx = 0.5d0*(pup + pdown)
          else
           xx = 0.5d0*(pup + pdown)
          end if
 
! Calculate Legendre polynomial from recursion relation and record the number
! of sign changes.
          pl1 = 1.0d0
          pl2 = xx
          node = 0
          do itheta = 2, ntheta_max
           pl3 = (2.0d0*real(itheta) - 1.0d0)*xx*pl2                         &
     &          - (real(itheta) - 1.0d0)*pl1
           if (sign(1.0d0,pl3) .ne. sign(1.0d0,pl2)) node = node + 1
           pl1 = pl2
           pl2 = pl3/real(itheta)
          end do

! Calculate the derivative
          plp = ntheta_max*(pl1 - xx*pl2)/(1.0d0 - xx*xx)
          if (node .ge. iterm) then
           pdown = xx
           ndown = node
          else
           pup = xx
           nup = node
          end if
          if (abs(pl2) .lt. 1.0d-10) exit 
         end do
         if (abs(pl2) .gt. 1.0d-10)                                          &
!     &    write (*,*) ' Warning no convergence in gleg after', jterm,        &
!     &                ' iterations'
         ctheta(iterm) = xx

! Gauss-Legendre weights
         ctheta_weights(iterm) = 2.0d0/((1.0d0 - xx*xx)*plp*plp)
        end do
 
! Fill in other half of points and weights from symmetry
        do iterm = 1, nw
         ctheta(ntheta_max - iterm + 1) = ctheta(iterm)
         ctheta(iterm) = - ctheta(iterm)
         ctheta_weights(ntheta_max - iterm + 1) = ctheta_weights(iterm)
        end do
        if (int(nw/2)*2 .ne. nw) ctheta(nw) = 0.0d0 ! Symmetry

! Format Statements
! ===========================================================================
        return
        end
 
