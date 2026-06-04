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
!                      Richard B. Evans
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
! Texas A&M - Traian Dumitrica
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


! density_calc.f
! Program Description
! ======================================================================
!
!  The subroutine density_calc will calculate the density and its
!  first and second derivatives for a given radius, and atom type
!
! ======================================================================
! Code written by:
! Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! Salt Lake City, UT 84112-0850
! (801) 581-8606 (office)      email: rbevans@hec.utah.edu
! (801) 581-4353 (fax)
!
! ======================================================================
        subroutine density_calc (iexc, ix, ispec, itype1, itype2, &
     &                           itype3, r, dr, dens, ddens, dddens)
        use precision, only: wp
        implicit none
        include '../parameters.inc'
        include '../exchange.inc'

! Argument Declaration and Description
! ======================================================================
! Input
        integer iexc
        integer ispec
        integer itype1
        integer itype2
        integer itype3
        integer ix

        real(kind=wp) r
        real(kind=wp) dr

! Output
        real(kind=wp) dens
        real(kind=wp) ddens
        real(kind=wp) dddens

! Local Parameters and Data Declaration
! ===========================================================================
        integer j1(7)
        data j1 /1, 1, 1, 2, 2, 3, 3/

        integer j2(7)
        data j2 /0, -1, 1, -1, 1, -1, 1/

! Local Variable Delcaration and Descrition
! ======================================================================
        integer issh
        integer itype
        integer j1ch
        integer j1at
        integer jssh

        integer in (3)

        real(kind=wp) densdr
        real(kind=wp) dens2dr
        real(kind=wp) dens_dr
        real(kind=wp) psiofr
        real(kind=wp) xinv4pi

        external psiofr

! Procedure
! ===========================================================================
! Initialize 1/4pi
        xinv4pi = 1.0d0/(4.0d0*3.141592653589793238462643D0)

        in(1) = itype1
        in(2) = itype2
        in(3) = itype3
        itype = in(ispec)

! Here the density is computed for r, r+dr and r-dr
        dens = 0.0d0
        do issh = 1, nsshxc(itype)
         dens = dens + xnocc(issh,itype)*psiofr(itype,issh,r)**2
        end do

! Here the density with respect to the charge correction term is calculated
! and added to dens, densdr and dens2dr.  The variable switch determines
! whether the correction is for the one, two or three center case.
        j1ch = j1(ix)
        if (j1ch .eq. ispec) then
         j1at = in(j1ch)
         jssh = iderorb(j1at)
         dens = dens + j2(ix)*dqorb(j1at)*psiofr(j1at,jssh,r)**2
        end if


! Only calculate the derivatives if doing GGA exchange-correlation.
! ***************************************************************************
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
     &      .or. iexc .eq. 9 .or. iexc .eq. 10) then

         densdr = 0.0d0
         dens_dr = 0.0d0
         do issh = 1, nsshxc(itype)
          densdr = &
     &     densdr + xnocc(issh,itype)*psiofr(itype,issh,r+dr)**2
          dens_dr = &
     &     dens_dr + xnocc(issh,itype)*psiofr(itype,issh,r-dr)**2
         end do

         if (j1ch .eq. ispec) then
          densdr = &
     &     densdr + j2(ix)*dqorb(j1at)*psiofr(j1at,jssh,r+dr)**2
          dens_dr = &
     &     dens_dr + j2(ix)*dqorb(j1at)*psiofr(j1at,jssh,r-dr)**2
         end if

! Here the first and second derivatives of the density is computed.
         if ((r - dr) .gt. 1.0d-5) then
          ddens = (densdr - dens_dr)/(2.0d0*dr)
          dddens = (densdr - 2.0d0*dens + dens_dr)/dr**2
         else

! At the endpoint do a forward difference. First, we need the point at r+2dr.
          dens2dr = 0.0d0
          do issh = 1, nsshxc(itype)
           dens2dr = &
     &      dens2dr + xnocc(issh,itype)*psiofr(itype,issh,r+2.0d0*dr)**2
          end do

          if (j1ch .eq. ispec) then
           dens2dr = dens2dr &
     &      + j2(ix)*dqorb(j1at)*psiofr(j1at,jssh,r+2.0d0*dr)**2
          end if

          ddens = (densdr - dens)/dr
          dddens = (dens2dr - 2.0d0*densdr + dens)/dr**2
         end if
        end if

! Convert to the correct units
        dens = dens*xinv4pi
        ddens = ddens*xinv4pi
        dddens = dddens*xinv4pi

! Format Statements
! ===========================================================================

        return
        end