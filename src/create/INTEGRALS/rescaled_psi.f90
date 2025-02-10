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

 
! rescaled_psi.f90
! Program Description
! ===========================================================================
!       This routine rescales the wavefunction according to the quantum
! number associated with angular components. Therefore, the input is a
! purely radial wavefunction, but the output is a radial wavefunction with
! an angular component.
!
! ===========================================================================
! Code written by:
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
        function rescaled_psi (l, m, rho, r, z, psi)
        use precision
        implicit none
        real(kind=long) rescaled_psi
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: m
 
        real(kind=long), intent (in) :: r
        real(kind=long), intent (in) :: rho
        real(kind=long), intent (in) :: z
        real(kind=long), intent (in) :: psi
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
!                       >>> THE MAGIC FORMULA <<<
 
! The i'th matrix element XY_Z - where X,Y = s, p, d, or f and Z = sigma, pi,
! delta, or phi - is:
!
!      faktor(i)*integral[X_Z(1)*Y_Z(2)*R1*R2*rho*drho*dz],
!
! where the 1, 2 mean the orbital is on atom 1,2. The faktor(i) are given
! below (e.g., faktor(5) is for the matrix element pp_pi and is 3/4)
! and the X_Z and Y_Z are:
!
!          -------------------------------------------------------
!          |                Magic formula table                  |
!          -------------------------------------------------------
!          |             s(sigma)  = 1                           |
!          |             p_sigma   = z/r                         |
!          |             p_pi      = rho/r                       |
!          |             d_sigma   = (2*z**2-rho**2)/r**2        |
!          |             d_pi      = rho*z/r**2                  |
!          |             d_delta   = rho**2/r**2                 |
!          |             f_sigma   = z*(2*z**2-3*rho**2)/r**3    |
!          |             f_pi      = rho*(4*z**2-rho**2)/r**3    |
!          |             f_delta   = z*rho**2/r**3               |
!          |             f_phi     = rho**3/r**3                 |
!          -------------------------------------------------------
! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
! For the short-range coulomb interactions make spherically symmetric
! s-states:
        if (l .eq. 0) then
         rescaled_psi = psi
         return
        end if

! Quick returns
        if (psi .eq. 0.0d0) then
         rescaled_psi = 0
         return
        end if
        if (r .le. 1.0d-5)then
          rescaled_psi = 0 
          return
        end if
 
! p-states:
        if (l .eq. 1) then
         if (abs(m) .eq. 1) rescaled_psi = psi*rho/r
         if (m .eq. 0) rescaled_psi = psi*z/r
         return
        end if
 
! d-states:
        if (l .eq. 2) then
         if (abs(m) .eq. 2) rescaled_psi = psi*rho**2/r**2
         if (abs(m) .eq. 1) rescaled_psi = psi*rho*z/r**2
         if (m .eq. 0) rescaled_psi = psi*(2.0d0*z**2 - rho**2)/r**2
         return
        end if
 
! f states:
        if (l .eq. 3) then
         if (abs(m) .eq. 3) rescaled_psi = psi*rho**3/r**3
         if (abs(m) .eq. 2) rescaled_psi = psi*rho**2*z/r**3
         if (abs(m) .eq. 1) rescaled_psi = psi*rho*(4.0d0*z**2 - rho**2)/r**3
         if (m .eq. 0) rescaled_psi = psi*z*(2.0d0*z**2 - 3.0d0*rho**2)/r**3
         return
        end if
 
! These terms for l > 4 are only used in the case of exact exchange.
        if (l .eq. 4) then
         if (abs(m) .eq. 4) rescaled_psi = psi*rho**4/r**4
         if (abs(m) .eq. 3) rescaled_psi = psi*z*rho**3/r**4
         if (abs(m) .eq. 2) rescaled_psi = psi*(6.0d0*z**2 - rho**2)*rho**2/r**4
         if (abs(m) .eq. 1)                                                   &
     &    rescaled_psi = psi*z*(4.0d0*z**2 - 3.0d0*rho**2)*rho/r**4
         if (m .eq. 0)                                                        &
     &    rescaled_psi = psi*(35.0d0*z**4/r**4                                &
     &                        + (3.0d0*rho**2 - 27.0d0*z**2)/r**2)
         return
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end