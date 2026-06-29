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

! potentials.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of the (non)-neutral
!       potentials
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module potentials
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use :: constants, only: abohr
  use :: math, only: math_interp_t, math_interp_new
  use :: utils, only: utils_open
  implicit none
  private

  public :: pot_atoms, pot_init, pot_end

  type :: pot_orbital_t
    private
    real(dp) :: rcut, etot
    type(math_interp_t) :: fr
  contains
    private
    procedure :: end => pot_orbital_end
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
    procedure, public :: get_energy => pot_atom_energy
    procedure, public :: get_vnn => pot_atom_vnn
    procedure :: end => pot_atom_end
  end type pot_atom_t

  type(pot_atom_t), allocatable :: pot_atoms(:)

contains

  subroutine pot_init()
    integer :: i, j, io, io2, ish, nspec, nshells, nz
    character(80) :: inname
    character(1000) :: fpathtemp
    character(1000), allocatable :: fpaths(:)

    io = utils_open("create.input", "r")
    read (io, *) nspec
    allocate (pot_atoms(nspec))
    do i = 1, nspec
      read (io, *) inname
      io2 = utils_open(inname, "r")
      read (io2, *)
      read (io2, *) nz
      do j = 1, 2
        read (io2, *)
      end do
      read (io2, *) fpathtemp
      read (io2, *) nshells
      allocate (fpaths(0:nshells))
      fpaths(0) = fpathtemp
      do ish = 1, nshells
        do j = 1, 4
          read (io2, *)
        end do
        read (io2, *) fpaths(ish)
      end do
      close (io2)
      pot_atoms(i) = pot_new_atom(nz, nshells, fpaths)
      deallocate (fpaths)
    end do
    close (io)
  end subroutine pot_init

  type(pot_atom_t) function pot_new_atom(nz, nshells, fpaths)
    integer, intent(in) :: nz, nshells
    character(1000), intent(in) :: fpaths(0:)
    integer :: i
    type(pot_orbital_t), allocatable :: pots(:)
    allocate (pots(0:nshells))
    do i = 0, nshells
      pots(i) = pot_new_orbital(i, fpaths(i))
    end do
    pot_new_atom = pot_atom_t(nz=nz, nshells=nshells, pots=pots)
    deallocate (pots)
  end function pot_new_atom

  type(pot_orbital_t) function pot_new_orbital(ish, fpath)
    integer, intent(in) :: ish
    character(1000), intent(in) :: fpath
    integer :: i, np, io
    real(dp) :: rcut, etot
    real(dp), allocatable :: r(:), vnn(:)
    type(math_interp_t) :: fr
    io = utils_open(fpath, "r")
    read (io, *)
    read (io, *)
    read (io, *) rcut
    read (io, *) np
    etot = 0.0_dp
    if (ish == 0) read (io, *) etot
    rcut = rcut*abohr
    allocate (r(np), vnn(np))
    do i = 1, np
      read (io, "(2d24.16)") r(i), vnn(i)
    end do
    close (io)
    fr = math_interp_new(r, vnn)
    deallocate (r, vnn)
    pot_new_orbital = pot_orbital_t(rcut=rcut, etot=etot, fr=fr)
  end function pot_new_orbital

  subroutine pot_orbital_end(this)
    class(pot_orbital_t), intent(inout) :: this
    call this%fr%end()
  end subroutine pot_orbital_end

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
    pot_atom_rcut = this%pots(ish)%rcut
  end function pot_atom_rcut

  pure real(dp) function pot_atom_energy(this)
    class(pot_atom_t), intent(in) :: this
    pot_atom_energy = this%pots(0)%etot
  end function pot_atom_energy

  pure real(dp) function pot_atom_vnn(this, ish, r, order)
    class(pot_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    real(dp), intent(in) :: r
    integer, intent(in), optional :: order
    integer :: o
    if ((r <= 0.0_dp) .or. (r >= this%pots(ish)%rcut)) then
      pot_atom_vnn = 0.0_dp
      return
    end if
    o = 0
    if (present(order)) o = order
    select case (o)
    case (0)
      pot_atom_vnn = this%pots(ish)%fr%f(r)
    case (1)
      pot_atom_vnn = this%pots(ish)%fr%df(r)
    case (2)
      pot_atom_vnn = this%pots(ish)%fr%ddf(r)
    case default
      pot_atom_vnn = 0.0_dp
    end select
  end function pot_atom_vnn

  subroutine pot_atom_end(this)
    class(pot_atom_t), intent(inout) :: this
    integer :: ish
    do ish = 0, this%nshells
      call this%pots(ish)%end()
    end do
    deallocate (this%pots)
  end subroutine pot_atom_end

  subroutine pot_end()
    integer :: nspec, ispec
    nspec = size(pot_atoms)
    do ispec = 1, nspec
      call pot_atoms(ispec)%end()
    end do
    deallocate (pot_atoms)
  end subroutine pot_end

end module potentials
