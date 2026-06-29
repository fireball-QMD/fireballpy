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

! pseudopotentials.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of the (non)-neutral
!       pseudopotentials
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module pseudopotentials
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use :: constants, only:abohr
  use :: math, only: math_interp_t, math_interp_new
  use :: utils, only: utils_open
  implicit none
  private

  public :: pp_atoms, pp_init, pp_end

  type :: pp_orbital_t
    private
    integer :: l
    real(dp) :: rcut, cl
    type(math_interp_t) :: fr
  contains
    private
    procedure :: end => pp_orbital_end
  end type pp_orbital_t

  type :: pp_atom_t
    private
    integer :: nz, nshells, iexc
    real(dp) :: rcut
    type(pp_orbital_t), allocatable :: pps(:)
  contains
    private
    procedure, public :: get_nz => pp_atom_nz
    procedure, public :: get_nshells => pp_atom_nshells
    procedure, public :: get_angular_momentum => pp_atom_angular_momentum
    procedure, public :: get_iexc => pp_atom_iexc
    procedure, public :: get_rcut => pp_atom_rcut
    procedure, public :: get_cl => pp_atom_cl
    procedure, public :: get_vpp => pp_atom_vpp
    procedure :: end => pp_atom_end
  end type pp_atom_t

  type(pp_atom_t), allocatable :: pp_atoms(:)

contains

  subroutine pp_init()
    integer :: i, j, io, io2, nspec, nshells, nz
    character(80) :: inname
    character(1000) :: fpath

    io = utils_open("create.input", "r")
    read (io, *) nspec
    allocate (pp_atoms(nspec))
    do i = 1, nspec
      read (io, *) inname
      io2 = utils_open(inname, "r")
      read (io2, *)
      read (io2, *) nz
      read (io2, *)
      read (io2, *) fpath
      read (io2, *)
      read (io2, *) nshells
      close (io2)
      pp_atoms(i) = pp_new_atom(nz, nshells, fpath)
      if (pp_atoms(i)%iexc /= pp_atoms(1)%iexc) stop
    end do
    close (io)
  end subroutine pp_init

  type(pp_atom_t) function pp_new_atom(nz, nshells, fpath)
    integer, intent(in) :: nz, nshells
    character(1000), intent(in) :: fpath
    integer :: i, ish, io, iexc, np, l
    real(dp) :: rcut, cl
    real(dp), allocatable :: r(:), vpp(:)
    type(pp_orbital_t), allocatable :: pps(:)
    allocate (pps(nshells))
    io = utils_open(fpath, "r")
    do i = 1, 14
      read (io, *)
    end do
    read (io, *) iexc
    do i = 1, 4
      read (io, *)
    end do
    read (io, *) rcut
    read (io, *) np
    do i = 1, np
      read (io, *)
    end do
    do ish = 1, nshells
      read (io, "(12x,i5)") np
      do i = 1, np
        read (io, *)
      end do
    end do
    do ish = 1, nshells
      read (io, "(3x,i1,8x,i5,4x,f14.7)") l, np, cl
      allocate (r(np), vpp(np))
      do i = 1, np
        read (io, *) r(i), vpp(i)
      end do
      pps(ish) = pp_new_orbital(l, rcut, cl, r, vpp)
      deallocate (r, vpp)
    end do
    close (io)
    pp_new_atom = pp_atom_t(nz=nz, nshells=nshells, iexc=iexc, rcut=rcut, pps=pps)
    deallocate (pps)
  end function pp_new_atom

  type(pp_orbital_t) function pp_new_orbital(l, rcut, cl, r, vpp)
    integer, intent(in) :: l
    real(dp), intent(in) :: rcut, cl
    real(dp), intent(in) :: r(:), vpp(:)
    type(math_interp_t) :: fr
    fr = math_interp_new(r, vpp)
    pp_new_orbital = pp_orbital_t(l=l, rcut=rcut, cl=cl, fr=fr)
  end function pp_new_orbital

  subroutine pp_orbital_end(this)
    class(pp_orbital_t), intent(inout) :: this
    call this%fr%end()
  end subroutine pp_orbital_end

  pure integer function pp_atom_nz(this)
    class(pp_atom_t), intent(in) :: this
    pp_atom_nz = this%nz
  end function pp_atom_nz

  pure integer function pp_atom_nshells(this)
    class(pp_atom_t), intent(in) :: this
    pp_atom_nshells = this%nshells
  end function pp_atom_nshells

  pure integer function pp_atom_angular_momentum(this, ish)
    class(pp_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    pp_atom_angular_momentum = this%pps(ish)%l
  end function pp_atom_angular_momentum

  pure integer function pp_atom_iexc(this)
    class(pp_atom_t), intent(in) :: this
    pp_atom_iexc = this%iexc
  end function pp_atom_iexc

  pure real(dp) function pp_atom_rcut(this, ish)
    class(pp_atom_t), intent(in) :: this
    integer, intent(in), optional :: ish
    if (present(ish)) then
      pp_atom_rcut = this%pps(ish)%rcut
    else
      pp_atom_rcut = this%rcut
    end if
  end function pp_atom_rcut

  pure real(dp) function pp_atom_cl(this, ish)
    class(pp_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    pp_atom_cl = this%pps(ish)%cl
  end function pp_atom_cl

  pure real(dp) function pp_atom_vpp(this, ish, r, order)
    class(pp_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    real(dp), intent(in) :: r
    integer, intent(in), optional :: order
    integer :: o
    if ((r <= 0.0_dp) .or. (r >= this%pps(ish)%rcut)) then
      pp_atom_vpp = 0.0_dp
      return
    end if
    o = 0
    if (present(order)) o = order
    select case (o)
    case (0)
      pp_atom_vpp = this%pps(ish)%fr%f(r)
    case (1)
      pp_atom_vpp = this%pps(ish)%fr%df(r)
    case (2)
      pp_atom_vpp = this%pps(ish)%fr%ddf(r)
    case default
      pp_atom_vpp = 0.0_dp
    end select
  end function pp_atom_vpp

  subroutine pp_atom_end(this)
    class(pp_atom_t), intent(inout) :: this
    integer :: ish
    do ish = 1, this%nshells
      call this%pps(ish)%end()
    end do
    deallocate (this%pps)
  end subroutine pp_atom_end

  subroutine pp_end()
    integer :: nspec, ispec
    nspec = size(pp_atoms)
    do ispec = 1, nspec
      call pp_atoms(ispec)%end()
    end do
    deallocate (pp_atoms)
  end subroutine pp_end

end module pseudopotentials
