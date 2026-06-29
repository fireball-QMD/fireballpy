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
    procedure, public :: f => math_interp_f
    procedure, public :: df => math_interp_df
    procedure, public :: ddf => math_interp_ddf
    procedure, public :: end => math_interp_end
  end type math_interp_t

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

  pure type(math_interp_t) function math_interp_new(x, y)
    real(kind=dp), intent(in) :: x(:), y(:)
    real(kind=dp), allocatable :: dx(:), dy(:), b(:), d(:), z(:), coefs(:, :)
    integer :: i, np, np1, np2
    np = size(x)
    np1 = np - 1
    np2 = np - 2
    allocate (coefs(4, np1), dx(np1), dy(np1), b(np), d(np), z(np))
    dx = x(2:np) - x(1:np1)
    dy = y(2:np) - y(1:np1)

    b(1) = 2.0_dp
    b(2) = 3.5_dp
    b(3:np1) = 3.75_dp
    b(np) = 1.75_dp
    d(1) = 0.0_dp
    d(2:np1) = 3.0_dp*(y(3:np) - y(1:np2))
    d(np) = 0.0_dp
    do i = 3, np
      d(i) = d(i) - 0.25_dp*d(i - 1)
    end do

    z(np) = d(np)/b(np)
    do i = np1, 1, -1
      z(i) = (d(i) - z(i + 1))/b(i)
    end do
    coefs(1, :) = y(1:np1)
    coefs(2, :) = z(1:np1)/dx
    coefs(3, :) = (3.0_dp*dy - 2*z(1:np1) - z(2:np))/dx**2
    coefs(4, :) = (-2.0_dp*dy + z(1:np1) + z(2:np))/dx**3
    math_interp_new = math_interp_t(np=np, x=x, y=y, coefs=coefs)
    deallocate (dx, dy, b, d, z, coefs)
  end function math_interp_new

  pure integer function math_interp_get_index(this, x)
    class(math_interp_t), intent(in) :: this
    real(kind=dp), intent(in) :: x
    integer :: i
    math_interp_get_index = -1
    do i = 1, (this%np - 1)
      if ((x >= this%x(i)) .and. (x < this%x(i + 1))) then
        math_interp_get_index = i
        return
      end if
    end do
  end function math_interp_get_index

  pure real(kind=dp) function math_interp_f(this, x)
    class(math_interp_t), intent(in) :: this
    real(kind=dp), intent(in) :: x
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_f = 0.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_f = this%coefs(1, i) + &
    &               dx*(this%coefs(2, i) + dx*(this%coefs(3, i) + dx*this%coefs(4, i)))
  end function math_interp_f

  pure real(kind=dp) function math_interp_df(this, x)
    class(math_interp_t), intent(in) :: this
    real(kind=dp), intent(in) :: x
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_df = 0.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_df = this%coefs(2, i) + dx*(2.0_dp*this%coefs(3, i) + dx*3.0_dp*this%coefs(4, i))
  end function math_interp_df

  pure real(kind=dp) function math_interp_ddf(this, x)
    class(math_interp_t), intent(in) :: this
    real(kind=dp), intent(in) :: x
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_ddf = 0.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_ddf = 2.0_dp*this%coefs(3, i) + dx*6.0_dp*this%coefs(4, i)
  end function math_interp_ddf

  subroutine math_interp_end(this)
    class(math_interp_t), intent(inout) :: this
    deallocate (this%x, this%y, this%coefs)
  end subroutine math_interp_end
end module math
