! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
! University of Utah - James P. Lewis, Chair
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
! Ohio State University - Dave Drabold

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

 
! gauntReal.f90
! Program Description
! ===========================================================================
!       This routine calculates:
!       sum_m Integral(Y*lm Xl2m2 Xl3m3)*Integral(Ylm Xl3m3 Xl4m4)
!       X the form of the real spherical harmonics
!         = i/sqrt(2)*(Ylm-(-1)**mY-ml) m < 0
!       X = Y0l                         m = 0
!         = 1/sqrt(2)*(Y-ml+(-1)**mYml) m > 0 
!       
!       and Gaunt coefficient  is defined as the integral over three
!       spherical harmonics 
!       gaunt(l1,l2,l3,m1,m2,m3)=Integral(Yl1m1 Yl2m2 Yl3m3)   
! ===========================================================================
! Code originally written by:
! Code written by:
! Daniel G. Trabada
! Departamento de físca teórica de la materia condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real(kind=wp) function gauntReal (l,l1,l2,l3,l4,m1,m2,m3,m4)
        use precision, only: wp
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: l1
        integer, intent (in) :: l2
        integer, intent (in) :: l3
        integer, intent (in) :: l4
        integer, intent (in) :: m1
        integer, intent (in) :: m2
        integer, intent (in) :: m3
        integer, intent (in) :: m4

! Local Parameters and Data Declaration
! ===========================================================================
        integer :: m
        complex(kind=wp) :: YXX
        complex(kind=wp) :: gauntComplex
 
! Procedure
! ===========================================================================
        gauntComplex=(0,0)
        do m=-l,l
          gauntComplex = gauntComplex + ((-1)**m)*YXX(l,l1,l2,-m,m1,m2)*YXX(l,l3,l4,m,m3,m4)
        end do
        if (abs(aimag(gauntComplex)) .gt. 1.0d-04 ) then 
         write (*,*) 'sum_m {YXX*YXX}, l1,l2,m1,m2,l3,l4,m3,m4 =',l1,l2,m1,m2,l3,l4,m3,m4 
         write (*,*) 'Warning function gauntReal it is not real',gauntComplex
        end if
        gauntReal = REAL(gauntComplex)
        return 
        end function gauntReal

! Program Declaration
! ===========================================================================
        complex(kind=wp) function YXX (l,l3,l4,m,m3,m4)
        use precision, only: wp
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: l3
        integer, intent (in) :: l4
        integer, intent (in) :: m
        integer, intent (in) :: m3
        integer, intent (in) :: m4
        real(kind=wp) :: aux
        real(kind=wp) :: gaunt
! Procedure
! ===========================================================================
        YXX=(0,0)
        if (m3 < 0 .and. m4 < 0) then 
          aux=-0.50d0*(gaunt(l,l3,l4,m,m3,m4)-((-1)**m4)*gaunt(l,l3,l4,m,m3,-m4) &
          &   -((-1)**m3)*gaunt(l,l3,l4,m,-m3,m4)+((-1)**(m3+m4))*gaunt(l,l3,l4,m,-m3,-m4))
          YXX = cmplx(aux,0.0d0)
        endif 
        if (m3 < 0 .and. m4 == 0) then
          aux=(sqrt(0.5d0))*(gaunt(l,l3,l4,m,m3,m4)-((-1)**m3)*gaunt(l,l3,l4,m,-m3,m4))
          YXX = cmplx(0.0d0,aux)
        endif 
        if (m3 < 0 .and. m4 > 0) then
          aux = 0.5*(gaunt(l,l3,l4,m,m3,-m4)+((-1)**m4)*gaunt(l,l3,l4,m,m3,m4)     &
          &   -((-1)**m3)*gaunt(l,l3,l4,m,-m3,-m4)-((-1)**(m3+m4))*gaunt(l,l3,l4,m,-m3,m4))
          YXX = cmplx(0.0d0,aux)
        endif 
        if (m3 == 0 .and. m4 < 0) then
          aux = (sqrt(0.5d0))*(gaunt(l,l3,l4,m,m3,m4)-((-1)**m4)*gaunt(l,l3,l4,m,m3,-m4))
          YXX = cmplx(0.0d0,aux)
        endif 
        if (m3 == 0 .and. m4 == 0) then
          aux = gaunt(l,l3,l4,m,m3,m4)
          YXX = cmplx(aux,0.0d0)
        endif 
         if (m3 == 0 .and. m4 > 0) then
          aux = (sqrt(0.5d0))*(gaunt(l,l3,l4,m,m3,-m4)+((-1)**m4)*gaunt(l,l3,l4,m,m3,m4))
          YXX = cmplx(aux,0.0d0)
         endif 
         if (m3 > 0 .and. m4 < 0) then
           aux = 0.5*(gaunt(l,l3,l4,m,-m3,m4)-((-1)**m4)*gaunt(l,l3,l4,m,-m3,-m4)     &
           &   +((-1)**m3)*gaunt(l,l3,l4,m,m3,m4)-((-1)**(m3+m4))*gaunt(l,l3,l4,m,m3,-m4))
           YXX = cmplx(0.0d0,aux)
         endif 
         if (m3 > 0 .and. m4 == 0) then
           aux = (sqrt(0.5d0))*(gaunt(l,l3,l4,m,-m3,m4)+(-1*m3)*gaunt(l,l3,l4,m,m3,m4))
           YXX = cmplx(aux,0.0d0)
         endif 
         if (m3 > 0 .and. m4 > 0) then
           aux = 0.5*(gaunt(l,l3,l4,m,-m3,-m4)+((-1)**m4)*gaunt(l,l3,l4,m,-m3,m4)      &
           & +((-1)**m3)*gaunt(l,l3,l4,m,m3,-m4)+((-1)**(m3+m4))*gaunt(l,l3,l4,m,m3,m4))
           YXX = cmplx(aux,0.0d0)
         endif 
         return 
      end function YXX
