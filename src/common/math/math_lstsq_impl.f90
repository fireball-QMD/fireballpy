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
submodule (math) math_lstsq_impl
implicit none

contains

module procedure dlstsq1
  real(kind=dp), allocatable :: newb(:, :)
  allocate (newb(size(b), 1))
  newb(:, 1) = b
  dlstsq1 = dlstsq2(a, newb, is_a_trans=is_a_trans)
  b = newb(:, 1)
end procedure dlstsq1

module procedure dlstsq2
  integer :: nrows, ncols, nrhs, nobs, lwork
  character(1) :: ctrans
  real(kind=dp), allocatable :: work(:)

  nrows = size(a, 1)
  ncols = size(a, 2)
  nobs = size(b, 1)
  nrhs = size(b, 2)
  ctrans = "N"
  if (present(is_a_trans)) then
    if (is_a_trans) ctrans = "T"
  end if
  if (((ctrans == "N") .and. (nrows /= nobs)) .or. ((ctrans == "T") .and. (ncols /= nobs))) then
    dlstsq2 = -1
    return
  end if
  allocate(work(1))
  call dgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, -1, dlstsq2)
  deallocate(work)
  allocate(work(lwork))
  call dgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, dlstsq2)
  if (dlstsq2 /= 0) return
  dlstsq2 = 0
end procedure dlstsq2

end submodule math_lstsq_impl
