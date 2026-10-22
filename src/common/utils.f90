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
!       This module exports utility functions for convenience
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module utils
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use :: constants, only: err_len
implicit none
private
public :: utils_open, utils_itos, utils_log, utils_nlines

integer, parameter, public :: UTILS_LOG_INFO = 0
integer, parameter, public :: UTILS_LOG_WARN = 1
integer, parameter, public :: UTILS_LOG_ERROR = 2

interface
  module integer function utils_open(filename, mode, skip)
    implicit none
    character(*), intent(in) :: filename
    character(*), intent(in), optional :: mode
    integer, intent(in), optional :: skip
  end function utils_open
end interface

interface
  module function utils_itos(i, lpad)
    implicit none
    integer, intent(in) :: i
    integer, intent(in), optional :: lpad
    character(:), allocatable, intent(out) :: utils_itos
  end function utils_itos
end interface

interface
  module subroutine utils_log(level, msg)
  integer, intent(in) :: level
  character(err_len), intent(in) :: msg
  end subroutine utils_log
end interface

interface
  module integer function utils_nlines(filename)
    implicit none
    character(*), intent(in) :: filename
  end function utils_nlines
end interface

end module utils
