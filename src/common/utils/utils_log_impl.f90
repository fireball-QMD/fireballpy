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
submodule (utils) utils_log_impl
implicit none

contains

module procedure utils_log
  select case (level)
  case (UTILS_LOG_INFO)
    write (stderr, "('[INFO] ',a)") trim(msg)
  case (UTILS_LOG_WARN)
    write (stderr, "('[WARN] ',a)") trim(msg)
  case (UTILS_LOG_ERROR)
    write (stderr, "('[ERROR] ',a)") trim(msg)
  case default
    write (stderr, "(a)") "[ERROR] incorrect level supplied to utils_log"
  end select
end procedure utils_log

end submodule utils_log_impl
