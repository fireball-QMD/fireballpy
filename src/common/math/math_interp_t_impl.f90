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

  module procedure math_interp_new
    integer :: i, np, np1, np2
    real(kind=dp), allocatable :: dx(:), dy(:), b(:), d(:), z(:), coefs(:, :)
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
    this%y = this%y * s
    this%coefs = this%coefs*s
  end procedure math_interp_rescale

  module procedure math_interp_end
    deallocate (this%x, this%y, this%coefs)
  end procedure math_interp_end
end submodule math_interp_t_impl
