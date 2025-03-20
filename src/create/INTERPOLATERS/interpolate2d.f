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


! interpolate2d.f
! Program Description
! ===========================================================================
!       This routine is a two-dimensional interpolater on a 4x4 sub-grid.
!
! ===========================================================================
! Code rewritten by:
! Kurt R. Glaesemann
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
        subroutine interpolate2d (rhoin, rhomin, rhomax, drho, zin,
     1                            zmin, zmax, dz, nnrho, nrho_points,
     2                            nnz, nz_points, frho, answer)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer nnrho
        integer nnz
        integer nrho_points
        integer nz_points

        real*8 drho
        real*8 dz
        real*8 rhoin
        real*8 rhomax
        real*8 rhomin
        real*8 zin
        real*8 zmax
        real*8 zmin

! This is the function being interpolated
        real*8 frho (nrho_points, nz_points)

! Output
        real*8 answer

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imidrho
        integer imidz
        integer j
        integer jj
        integer k

        real*8 e36t
        real*8 f0p3
        real*8 f0p6
        real*8 f1m2
        real*8 f1m3
        real*8 f1p3
        real*8 f1p6
        real*8 ftp
        real*8 gradrho
        real*8 gradtest
        real*8 gradz
        real*8 prho
        real*8 prod
        real*8 pz
        real*8 tp

        real*8 b (0:5)
        real*8 bb (0:5, -2:3)
        real*8 g (-2:3)
        real*8 fun (-1:2, -1:2)

! Procedure
! ===========================================================================
! Check and make sure that the point (r, z) is within the limits of the
! stored density.
        if (rhoin .gt. (rhomax + 1.0d-5) .or. rhoin .lt. rhomin) then
         write (*,*) ' What the heck is going on in interpolate2d.f !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input rhoin = ', rhoin
         write (*,*) ' Max. data = ', rhomax, ' Min. data = ', rhomin
        end if

        if (zin .gt. (zmax + 1.0d-5) .or. zin .lt. (zmin - 1.0d-5)) then
         write (*,*) ' What the heck is going on in interpolate2d.f !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input zin = ', zin
         write (*,*) ' Max. data = ', zmax, ' Min. data = ', zmin
        end if

! Much of the code in this routine is a dodge to avoid interpolating, even
! though the mission of this subprogram is to interpolate. The dodge can be
! approached at two levels: first, if the gradient is locally very small, we
! can just return. This is CRUCIAL for multirc problems, but seems to help for
! monorc.

! Set-up everything for interpolation.
        imidrho = int((rhoin - rhomin)/drho) + 1
        imidz = int((zin - zmin)/dz) + 1

        if (imidrho .lt. 2) imidrho = 2
        if (imidrho .gt. nnrho) imidrho = nnrho
        if (imidz .lt. 2) imidz = 2
        if (imidz .gt. nnz) imidz = nnz

        prho = (rhoin - rhomin)/drho - dfloat(imidrho - 1)
        pz = (zin - zmin)/dz - dfloat(imidz - 1)

        fun(-1,-1) = frho(imidrho - 1,imidz - 1)
        fun(-1, 0) = frho(imidrho - 1,imidz)
        fun(-1, 1) = frho(imidrho - 1,imidz + 1)
        fun(-1, 2) = frho(imidrho - 1,imidz + 2)

        fun(0,-1) = frho(imidrho,imidz - 1)
        fun(0, 0) = frho(imidrho,imidz)
        fun(0, 1) = frho(imidrho,imidz + 1)
        fun(0, 2) = frho(imidrho,imidz + 2)

        fun(1,-1) = frho(imidrho + 1,imidz - 1)
        fun(1, 0) = frho(imidrho + 1,imidz)
        fun(1, 1) = frho(imidrho + 1,imidz + 1)
        fun(1, 2) = frho(imidrho + 1,imidz + 2)

        fun(2,-1) = frho(imidrho + 2,imidz - 1)
        fun(2, 0) = frho(imidrho + 2,imidz)
        fun(2, 1) = frho(imidrho + 2,imidz + 1)
        fun(2, 2) = frho(imidrho + 2,imidz + 2)

! **************************************************************************
! If the gradient is small, then do quadratic quick bivariate interpolation.
        gradrho = (fun(0,0) - fun(1,0))/drho
        gradz = (fun(0,0) - fun(0,1))/dz

        gradtest = abs(gradrho) + abs(gradz)

! Form the criterion for a quick interpolation. Empirically, I find that
! gradtest < 1.0d-05 is adequate. If you dont want this option, change
! 1.0d-05 in the next line to 0.0d0!
        if (gradtest .lt. 1.0d-05) then
         answer = (1.0d0 - prho - pz)*fun(0,0) + prho*fun(1,0)
     1                                         + pz*fun(0,1)
         return
        end if

! **************************************************************************
! Phase III. All else fails. Interpolate carefully. Original pfed
! interpolator with minimal multiplies.
        e36t = 1.0d0/36.0d0

        do k = - 1, 2
         f1m2 = fun(k,-1) + fun(k,-1)
         f1m3 = f1m2 + fun(k,-1)

         f0p3 = fun(k,0) + fun(k,0) + fun(k,0)
         f0p6 = f0p3 + f0p3

         f1p3 = fun(k,1) + fun(k,1) + fun(k,1)
         f1p6 = f1p3 + f1p3

         bb(3,k) = - fun(k,-1) + f0p3 - f1p3 + fun(k,2)
         bb(2,k) = f1m3 - f0p6 + f1p3
         bb(1,k) = - f1m2 - f0p3 + f1p6 - fun(k,2)

         tp = fun(k,0)
         tp = tp + tp + tp
         bb(0,k) = tp + tp

         prod = bb(3,k)*pz
         do j = 1, 2
          jj = 3 - j
          ftp = bb(jj,k)
          prod = (prod + ftp)*pz
         end do
         g(k) = prod + bb(0,k)
        end do

        f1m2 = g(-1) + g(-1)
        f1m3 = f1m2 + g(-1)

        f0p3 = g(0) + g(0) + g(0)
        f0p6 = f0p3 + f0p3

        f1p3 = g(1) + g(1) + g(1)
        f1p6 = f1p3 + f1p3

        b(3) = -g(-1) + f0p3 - f1p3 + g(2)
        b(2) = f1m3 - f0p6 + f1p3
        b(1) = -f1m2 - f0p3 + f1p6 - g(2)
        tp = g(0) + g(0) + g(0)
        b(0) = tp + tp

        prod = b(3)*prho
        do j = 1, 2
         jj = 3 - j
         ftp = b(jj)
         prod = (prod + ftp)*prho
        end do
        prod = prod + b(0)

! Final answer
        answer = e36t*prod

! Format Statements
! ===========================================================================
 
        return
        end
