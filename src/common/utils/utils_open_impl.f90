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

! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
submodule (utils) utils_open_impl
implicit none

contains

module procedure utils_open
  integer :: i, stat
  character(3) :: mode_
  character(err_len) :: errmsg
  character(:), allocatable :: action_, position_, status_, access_, form_

  if (.not. present(mode)) then
    mode_ = "r t"
  else
    select case (len(mode))
    case (1)
      mode_ = mode//" t"
    case (2)
      mode_ = mode//"t"
    case default
      mode_ = trim(adjustl(mode(1:3)))
    end select
  end if

  select case (mode_(1:2))
  case ("r")
    action_ = "read"
    position_ = "asis"
    status_ = "old"
  case ("w")
    action_ = "write"
    position_ = "asis"
    status_ = "replace"
  case ("a")
    action_ = "write"
    position_ = "append"
    status_ = "old"
  case ("x")
    action_ = "write"
    position_ = "asis"
    status_ = "new"
  case ("r+")
    action_ = "readwrite"
    position_ = "asis"
    status_ = "old"
  case ("w+")
    action_ = "readwrite"
    position_ = "asis"
    status_ = "replace"
  case ("a+")
    action_ = "readwrite"
    position_ = "append"
    status_ = "old"
  case ("x+")
    action_ = "readwrite"
    position_ = "asis"
    status_ = "new"
  case default
    call utils_log(UTILS_LOG_ERROR, "Unknown mode '"//mode_(1:2)//"-'")
    error stop 1
  end select

  select case (mode_(3:3))
  case ("t")
    form_ = "formatted"
    access_ = "sequential"
  case ("b")
    form_ = "unformatted"
    access_ = "stream"
  case default
    call utils_log(UTILS_LOG_ERROR, "Unknown mode '--"//mode_(3:3)//"'")
    error stop 1
  end select

  open (newunit=utils_open, file=trim(filename), action=action_, position=position_, status=status_, &
      access=access_, form=form_, iostat=stat, iomsg=errmsg)
  if (stat /= 0) then
    call utils_log(UTILS_LOG_ERROR, achar(iachar(errmsg(1:1))+32)//trim(errmsg(2:)))
    error stop 1
  end if

  if (present(skip)) then
    do i = 1, skip
      read (utils_open, *)
    end do
  end if
end procedure utils_open

end submodule utils_open_impl
