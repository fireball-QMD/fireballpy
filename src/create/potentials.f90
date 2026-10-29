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
! adp with this program.  If not, see <http://www.gnu.org/licenses/>.

! potentials.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of the (non-)neutral
!       potentials
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module potentials
use, intrinsic :: iso_fortran_env, only: dp => real64
use, intrinsic :: ieee_arithmetic, only: ieee_signaling_nan, ieee_value
use :: constants, only: abohr, tolerance, PATH_LEN, LNAMES
use :: math, only: math_interp_t, math_interp_new
use :: utils, only: utils_open, utils_itos
use :: input, only: input_t
implicit none
private

public :: pot_atoms, pot_init, pot_end

type, extends(math_interp_t) :: pot_orbital_t
  private
  real(dp) :: rcut
end type pot_orbital_t

type :: pot_atom_t
  private
  integer :: nz, nshells
  type(pot_orbital_t), allocatable :: pots(:)
contains
  private
  procedure, public :: get_nz => pot_atom_nz
  procedure, public :: get_nshells => pot_atom_nshells
  procedure, public :: get_rcut => pot_atom_rcut
  procedure, public :: get_vneutral => pot_atom_vneutral
  procedure, public :: get_vnn => pot_atom_vnn
  procedure :: end => pot_atom_end
end type pot_atom_t

type(pot_atom_t), allocatable :: pot_atoms(:)

contains

subroutine pot_init(inp)
  type(input_t), intent(in) :: inp
  integer :: i, ish, nsh
  character(PATH_LEN), allocatable :: fpaths(:)
  allocate(pot_atoms(size(inp%species)))
  do i = 1, size(inp%species)
    associate (sp => inp%species(i))
      nsh = sp%nshells
      if (allocated(fpaths)) deallocate(fpaths)
      allocate(fpaths(0:nsh))
      fpaths(0) = trim(inp%folder)//trim(utils_itos(sp%znuc, lpad=3))//".na0"
      do ish = 1, nsh
        if (sp%shells(ish)%excited == 0) then
          fpaths(ish) = trim(inp%folder)//trim(utils_itos(sp%znuc, lpad=3))//"_"//LNAMES(sp%shells(ish)%l)//".na"
        else
          fpaths(ish) = trim(inp%folder)//trim(utils_itos(sp%znuc, lpad=3))//"_"//LNAMES(sp%shells(ish)%l)//".ena"
        end if
      end do
      pot_atoms(i) = pot_new_atom(sp%znuc, nsh, fpaths)
    end associate
  end do
end subroutine pot_init

type(pot_atom_t) function pot_new_atom(nz, nsh, fpaths)
  integer, intent(in) :: nz, nsh
  character(PATH_LEN), intent(in) :: fpaths(0:)
  integer :: i
  type(pot_orbital_t), allocatable :: pots(:)
  allocate(pots(0:nsh))
  do i = 0, nsh
    pots(i) = pot_new_orbital(i, fpaths(i))
  end do
  pot_new_atom = pot_atom_t(nz=nz, nshells=nsh, pots=pots)
end function pot_new_atom

type(pot_orbital_t) function pot_new_orbital(ish, fpath)
  integer, intent(in) :: ish
  character(1000), intent(in) :: fpath
  integer :: i, np, io
  real(dp) :: rcut
  real(dp), allocatable :: r(:), vnn(:)
  type(math_interp_t) :: fr
  io = utils_open(fpath, "r", skip=1)
  read (io, "(6x,i8)") np
  allocate(r(np), vnn(np))
  do i = 1, np
    read (io, "(2ES16.8)") r(i), vnn(i)
  end do
  close (io)
  fr = math_interp_new(r, vnn, rval=0.0_dp)
  rcut = r(np)
  pot_new_orbital = pot_orbital_t(rcut=rcut, math_interp_t=fr)
end function pot_new_orbital

pure integer function pot_atom_nz(this)
  class(pot_atom_t), intent(in) :: this
  pot_atom_nz = this%nz
end function pot_atom_nz

pure integer function pot_atom_nshells(this)
  class(pot_atom_t), intent(in) :: this
  pot_atom_nshells = this%nshells
end function pot_atom_nshells

pure real(dp) function pot_atom_rcut(this, ish)
  class(pot_atom_t), intent(in) :: this
  integer, intent(in) :: ish
  if (ish < 1 .or. ish > this%nshells) then
    pot_atom_rcut = ieee_value(pot_atom_rcut, ieee_signaling_nan)
    return
  end if
  pot_atom_rcut = this%pots(ish)%rcut
end function pot_atom_rcut

pure real(dp) function pot_atom_vneutral(this, r, order)
  class(pot_atom_t), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in), optional :: order
  pot_atom_vneutral = this%pots(0)%f(r, order=order)
end function pot_atom_vneutral

pure real(dp) function pot_atom_vnn(this, ish, r, order)
  class(pot_atom_t), intent(in) :: this
  integer, intent(in) :: ish
  real(dp), intent(in) :: r
  integer, intent(in), optional :: order
  integer :: o
  o = 0
  if (present(order)) o = order
  if (r >= this%pots(ish)%rcut) then
    select case (o)
    case (0)
      pot_atom_vnn = 1.0_dp/r
    case (1)
      pot_atom_vnn = -1.0_dp/(r*r)
    case (2)
      pot_atom_vnn = 2.0_dp/(r*r*r)
    case default
      pot_atom_vnn = ieee_value(pot_atom_vnn, ieee_signaling_nan)
    end select
    return
  end if
  pot_atom_vnn = this%pots(ish)%f(r, order=o)
end function pot_atom_vnn

subroutine pot_atom_end(this)
  class(pot_atom_t), intent(inout) :: this
  integer :: ish
  do ish = 0, this%nshells
    call this%pots(ish)%end()
  end do
  deallocate(this%pots)
end subroutine pot_atom_end

subroutine pot_end()
  integer :: nspec, ispec
  nspec = size(pot_atoms)
  do ispec = 1, nspec
    call pot_atoms(ispec)%end()
  end do
  deallocate(pot_atoms)
end subroutine pot_end

end module potentials
