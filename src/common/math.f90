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
! Arizona State University - Gary B. Adams Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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
! awp with this program.  If not, see <http://www.gnu.org/licenses/>.

! math.f90
! Program Description
! ==============================================================================
!       This module calculates handles mathematical algorithms
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module math
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: math_lstsq
  public :: math_interp_t, math_interp_new

  type :: math_interp_t
    private
    integer :: np
    real(dp), allocatable :: x(:), y(:), coefs(:, :)
  contains
    private
    procedure :: get_index => math_interp_get_index
    procedure, public :: get_n => math_interp_get_n
    procedure, public :: get_x => math_interp_get_x
    procedure, public :: get_y => math_interp_get_y
    procedure, public :: get_coef => math_interp_get_coef
    procedure, public :: f => math_interp_f
    procedure, public :: df => math_interp_df
    procedure, public :: ddf => math_interp_ddf
    procedure, public :: rescale => math_interp_rescale
    procedure, public :: end => math_interp_end
  end type math_interp_t

  interface
    module pure type(math_interp_t) function math_interp_new(x, y)
      implicit none
      real(kind=dp), intent(in) :: x(:), y(:)
      real(kind=dp), allocatable :: dx(:), dy(:), b(:), d(:), z(:), coefs(:, :)
    end function math_interp_new

    module pure integer function math_interp_get_index(this, x)
      implicit none
      class(math_interp_t), intent(in) :: this
      real(kind=dp), intent(in) :: x
    end function math_interp_get_index

    module pure integer function math_interp_get_n(this)
      implicit none
      class(math_interp_t), intent(in) :: this
    end function math_interp_get_n

    module pure real(kind=dp) function math_interp_get_x(this, i)
      implicit none
      class(math_interp_t), intent(in) :: this
      integer, intent(in) :: i
    end function math_interp_get_x

    module pure real(kind=dp) function math_interp_get_y(this, i)
      implicit none
      class(math_interp_t), intent(in) :: this
      integer, intent(in) :: i
    end function math_interp_get_y

    module pure real(kind=dp) function math_interp_get_coef(this, i, j)
      implicit none
      class(math_interp_t), intent(in) :: this
      integer, intent(in) :: i, j
    end function math_interp_get_coef

    module pure real(kind=dp) function math_interp_f(this, x)
      implicit none
      class(math_interp_t), intent(in) :: this
      real(kind=dp), intent(in) :: x
    end function math_interp_f

    module pure real(kind=dp) function math_interp_df(this, x)
      implicit none
      class(math_interp_t), intent(in) :: this
      real(kind=dp), intent(in) :: x
    end function math_interp_df

    module pure real(kind=dp) function math_interp_ddf(this, x)
      implicit none
      class(math_interp_t), intent(in) :: this
      real(kind=dp), intent(in) :: x
    end function math_interp_ddf

    module subroutine math_interp_rescale(this, s)
      implicit none
      class(math_interp_t), intent(inout) :: this
      real(kind=dp), intent(in) :: s
    end subroutine math_interp_rescale

    module subroutine math_interp_end(this)
      implicit none
      class(math_interp_t), intent(inout) :: this
    end subroutine math_interp_end

  end interface

  interface math_lstsq
    procedure :: lstsq1
    procedure :: lstsq2
  end interface math_lstsq

contains

  integer function lstsq1(a, b, is_a_trans)
    real(kind=dp), intent(inout) :: a(:, :), b(:)
    logical, intent(in), optional :: is_a_trans
    real(kind=dp), allocatable :: newb(:, :)
    allocate (newb(size(b), 1))
    newb(:, 1) = b
    if (present(is_a_trans)) then
      lstsq1 = lstsq2(a, newb, is_a_trans=is_a_trans)
    else
      lstsq1 = lstsq2(a, newb)
    end if
    b = newb(:, 1)
    deallocate (newb)
  end function lstsq1

  integer function lstsq2(a, b, is_a_trans)
    real(kind=dp), intent(inout) :: a(:, :), b(:, :)
    logical, intent(in), optional :: is_a_trans
    integer :: nrows, ncols, nrhs, nobs, lwork
    character(1) :: ctrans
    real(kind=dp), allocatable :: work(:)

    nrows = size(a, 1)
    ncols = size(a, 2)
    nobs = size(b, 1)
    nrhs = size(b, 2)
    if (.not. present(is_a_trans)) then
      ctrans = "N"
    else if (is_a_trans) then
      ctrans = "T"
    else
      ctrans = "N"
    end if
    if (((ctrans == "N") .and. (nrows /= nobs)) .or. ((ctrans == "T") .and. (ncols /= nobs))) then
      lstsq2 = -1
      return
    end if

    lwork = min(nrows, ncols) + max(min(nrows, ncols), nrhs)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate (work(lwork))
    call dgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, lstsq2)
    deallocate (work)
    if (lstsq2 /= 0) return
    lstsq2 = 0
  end function lstsq2
end module math
