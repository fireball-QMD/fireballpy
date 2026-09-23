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

! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
submodule (math) math_interp_t_impl
  implicit none

contains

  subroutine math_solve_tridiag(x, upper, lower, diag)
    real(kind=dp), intent(inout) :: x(:)
    real(kind=dp), intent(in) :: upper(:), lower(:), diag(:)
    integer :: i, i1, np, np1
    real(kind=dp) :: temp
    real(kind=dp), allocatable :: scratch(:)
    np = size(x)
    np1 = np - 1
    allocate (scratch(np1))

    scratch(1) = upper(1)/diag(1)
    x(1) = x(1)/diag(1)
    do i = 2, np1
      i1 = i - 1
      temp = 1.0_dp/(diag(i) - lower(i1)*scratch(i1))
      scratch(i) = upper(i)*temp
      x(i) = (x(i) - lower(i1)*x(i1))*temp
    end do
    x(np) = (x(np) - lower(np1)*x(np1))/(diag(np) - lower(np1)*scratch(np1))
    do i = np1, 1, -1
      x(i) = x(i) - scratch(i)*x(i + 1)
    end do
  end subroutine math_solve_tridiag

  module procedure math_interp_new
    integer :: i, np, np1, np2
    real(kind=dp), allocatable :: dx(:), s(:), a(:), b(:), c(:), d(:), temp(:), coefs(:, :)
    np = size(x)
    np1 = np - 1
    np2 = np - 2
    allocate (coefs(4, np1), dx(np1), s(np1), a(np1), b(np), c(np1), d(np), temp(np1))
    dx = x(2:) - x(:np1)
    s = (y(2:) - y(:np1))/dx

    ! Set tridiagonal coefs for middle equations
    a(:np2) = dx(2:)
    b(2:np1) = 2.0_dp*(dx(2:) + dx(:np2))
    c(2:) = dx(:np2)
    d(2:np1) = 3.0_dp*(dx(2:)*s(:np2) + dx(:np2)*s(2:))

    ! Left not-a-knot
    c(1) = x(3) - x(1)
    b(1) = dx(2)
    d(1) = ((dx(1) + 2.0_dp*c(1))*dx(2)*s(1) + dx(1)*dx(1)*s(2))/c(1)

    ! Right not-a-knot
    a(np1) = x(np) - x(np2)
    b(np) = dx(np2)
    d(np) = (dx(np1)*dx(np1)*s(np2) + (2.0_dp*a(np1) + dx(np1))*dx(np2)*s(np1))/a(np1)

    call math_solve_tridiag(d, c, a, b)

    temp = (d(:np1) + d(2:) - 2*s)/dx
    coefs(1, :) = y(:np1)
    coefs(2, :) = d(:np1)
    coefs(3, :) = (s - d(:np1))/dx - temp
    coefs(4, :) = temp/dx
    math_interp_new = math_interp_t(np=np, x=x, y=y, coefs=coefs)
  end procedure math_interp_new

  module procedure math_interp_get_index
    integer :: i
    math_interp_get_index = -1
    if (x < this%x(1)) return
    if (x >= this%x(this%np)) return
    do i = 1, (this%np - 1)
      if ((x >= this%x(i)) .and. (x < this%x(i + 1))) then
        math_interp_get_index = i
        return
      end if
    end do
    math_interp_get_index = -2
  end procedure math_interp_get_index

  module procedure math_interp_get_n
    math_interp_get_n = this%np
  end procedure math_interp_get_n

  module procedure math_interp_get_x
    math_interp_get_x = this%x(i)
  end procedure math_interp_get_x

  module procedure math_interp_get_y
    math_interp_get_y = this%y(i)
  end procedure math_interp_get_y

  module procedure math_interp_get_coef
    math_interp_get_coef = this%coefs(i, j)
  end procedure math_interp_get_coef

  module procedure math_interp_f
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_f = 0.0_dp
      return
    end if
    if (i == -2) then
      math_interp_f = 10000000000.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_f = this%coefs(1, i) + &
    &               dx*(this%coefs(2, i) + dx*(this%coefs(3, i) + dx*this%coefs(4, i)))
  end procedure math_interp_f

  module procedure math_interp_df
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_df = 0.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_df = this%coefs(2, i) + dx*(2.0_dp*this%coefs(3, i) + dx*3.0_dp*this%coefs(4, i))
  end procedure math_interp_df

  module procedure math_interp_ddf
    integer :: i
    real(kind=dp) :: dx
    i = this%get_index(x)
    if (i == -1) then
      math_interp_ddf = 0.0_dp
      return
    end if
    dx = x - this%x(i)
    math_interp_ddf = 2.0_dp*this%coefs(3, i) + dx*6.0_dp*this%coefs(4, i)
  end procedure math_interp_ddf

  module procedure math_interp_rescale
    this%y = s*this%y
    this%coefs = s*this%coefs
  end procedure math_interp_rescale

  module procedure math_interp_end
    deallocate (this%x, this%y, this%coefs)
  end procedure math_interp_end
end submodule math_interp_t_impl
