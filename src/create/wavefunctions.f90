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

! wavefunctions.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of the wavefunctions
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module wavefunctions
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use :: constants, only:abohr, inv4pi, tolerance
  use :: math, only: math_interp_t, math_interp_new
  use :: utils, only: utils_open
  implicit none
  private

  public :: wf_atoms, wf_dens, wf_init, wf_end

  interface wf_dens
    procedure :: dens1c
    procedure :: dens2c
  end interface wf_dens

  type :: wf_orbital_t
    private
    integer :: l
    real(dp) :: rcut, rcut_max, q, qref
    type(math_interp_t) :: fr
  contains
    private
    procedure :: end => wf_orbital_end
  end type wf_orbital_t

  type :: wf_atom_t
    private
    integer :: nz, nshells
    real(dp) :: rcut
    type(wf_orbital_t), allocatable :: wfs(:)
  contains
    private
    procedure, public :: get_nz => wf_atom_nz
    procedure, public :: get_nshells => wf_atom_nshells
    procedure, public :: get_angular_momentum => wf_atom_angular_momentum
    procedure, public :: get_charge => wf_atom_charge
    procedure, public :: get_ref_charge => wf_atom_ref_charge
    procedure, public :: get_rcut => wf_atom_rcut
    procedure, public :: get_psi => wf_atom_psi
    procedure :: end => wf_atom_end
  end type wf_atom_t

  type(wf_atom_t), allocatable :: wf_atoms(:)

contains

  subroutine wf_init()
    integer :: i, j, io, io2, ish, nspec, nshells, nz
    character(80) :: inname
    character(1000), allocatable :: fpaths(:)

    io = utils_open("create.input", "r")
    read (io, *) nspec
    allocate (wf_atoms(nspec))
    do i = 1, nspec
      read (io, *) inname
      io2 = utils_open(inname, "r")
      read (io2, *)
      read (io2, *) nz
      do j = 1, 3
        read (io2, *)
      end do
      read (io2, *) nshells
      allocate (fpaths(nshells))
      do ish = 1, nshells
        do j = 1, 3
          read (io2, *)
        end do
        read (io2, *) fpaths(ish)
        read (io2, *)
      end do
      close (io2)
      wf_atoms(i) = wf_new_atom(nz, nshells, fpaths)
      deallocate (fpaths)
    end do
    close (io)
  end subroutine wf_init

  type(wf_atom_t) function wf_new_atom(nz, nshells, fpaths)
    integer, intent(in) :: nz, nshells
    character(1000), intent(in) :: fpaths(:)
    integer :: i
    real(dp) :: rcut
    type(wf_orbital_t), allocatable :: wfs(:)
    allocate (wfs(nshells))
    rcut = 0.0_dp
    do i = 1, nshells
      wfs(i) = wf_new_orbital(fpaths(i))
      rcut = max(rcut, wfs(i)%rcut)
    end do
    wf_new_atom = wf_atom_t(nz=nz, nshells=nshells, rcut=rcut, wfs=wfs)
    deallocate (wfs)
  end function wf_new_atom

  type(wf_orbital_t) function wf_new_orbital(fpath)
    character(1000), intent(in) :: fpath
    integer :: i, j, np, remitems, nlines, l, io
    real(dp) :: dr, rcut, rcut_max, q
    real(dp), allocatable :: r(:), psi(:)
    type(math_interp_t) :: fr
    io = utils_open(fpath, "r")
    read (io, *)
    read (io, *)
    read (io, *) np
    read (io, *) rcut, rcut_max, q
    rcut = rcut*abohr
    read (io, *) l
    nlines = np/4
    remitems = mod(np, 4)
    allocate (r(np), psi(np))
    do i = 1, nlines
      read (io, "(4d18.10)") (psi(4*(i - 1) + j), j=1, 4)
    end do
    if (remitems /= 0) then
      read (io, "(4d18.10)") (psi(4*nlines + j), j=1, remitems)
    end if
    close (io)
    dr = rcut/real(np - 1, kind=dp)
    do i = 1, np
      r(i) = real(i - 1, kind=dp)*dr
    end do
    fr = math_interp_new(r, psi)
    deallocate (r, psi)
    ! TODO: allow qref to be different
    wf_new_orbital = wf_orbital_t(l=l, rcut=rcut, rcut_max=rcut_max, q=q, qref=q, fr=fr)
  end function wf_new_orbital

  subroutine wf_orbital_end(this)
    class(wf_orbital_t), intent(inout) :: this
    call this%fr%end()
  end subroutine wf_orbital_end

  pure integer function wf_atom_nz(this)
    class(wf_atom_t), intent(in) :: this
    wf_atom_nz = this%nz
  end function wf_atom_nz

  pure integer function wf_atom_nshells(this)
    class(wf_atom_t), intent(in) :: this
    wf_atom_nshells = this%nshells
  end function wf_atom_nshells

  pure integer function wf_atom_angular_momentum(this, ish)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    wf_atom_angular_momentum = this%wfs(ish)%l
  end function wf_atom_angular_momentum

  pure real(dp) function wf_atom_charge(this, ish)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    wf_atom_charge = this%wfs(ish)%q
  end function wf_atom_charge

  pure real(dp) function wf_atom_ref_charge(this, ish)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    wf_atom_ref_charge = this%wfs(ish)%qref
  end function wf_atom_ref_charge

  pure real(dp) function wf_atom_rcut(this, ish)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in), optional :: ish
    if (present(ish)) then
      wf_atom_rcut = this%wfs(ish)%rcut
    else
      wf_atom_rcut = this%rcut
    end if
  end function wf_atom_rcut

  pure real(dp) function wf_atom_psi(this, ish, r, order)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    real(dp), intent(in) :: r
    integer, intent(in), optional :: order
    integer :: o
    if ((r <= 0.0_dp) .or. (r >= this%wfs(ish)%rcut)) then
      wf_atom_psi = 0.0_dp
      return
    end if
    o = 0
    if (present(order)) o = order
    select case (o)
    case (0)
      wf_atom_psi = this%wfs(ish)%fr%f(r)
    case (1)
      wf_atom_psi = this%wfs(ish)%fr%df(r)
    case (2)
      wf_atom_psi = this%wfs(ish)%fr%ddf(r)
    case default
      wf_atom_psi = 0.0_dp
    end select
  end function wf_atom_psi

  subroutine wf_atom_end(this)
    class(wf_atom_t), intent(inout) :: this
    integer :: ish
    do ish = 1, this%nshells
      call this%wfs(ish)%end()
    end do
    deallocate (this%wfs)
  end subroutine wf_atom_end

  subroutine wf_end()
    integer :: nspec, ispec
    nspec = size(wf_atoms)
    do ispec = 1, nspec
      call wf_atoms(ispec)%end()
    end do
    deallocate (wf_atoms)
  end subroutine wf_end

  pure subroutine dens1c(ispec, dens, ddens, dddens)
    integer, intent(in) :: ispec
    real(dp), intent(out) :: dens(:)
    real(dp), intent(out), optional :: ddens(:, :), dddens(:, :, :)
    logical :: isder
    integer :: i, j, nr
    real(dp) :: r, dr, tpsi
    real(dp), allocatable :: psi(:)
    nr = size(dens)
    dr = wf_atoms(ispec)%get_rcut()/real(nr - 1, kind=dp)
    isder = present(ddens) .and. present(dddens)
    allocate (psi(wf_atoms(ispec)%get_nshells()))
    do i = 1, nr
      r = real(i - 1, kind=dp)*dr
      dens(i) = 0.0_dp
      if (r > tolerance) then
        do j = 1, wf_atoms(ispec)%get_nshells()
          psi(j) = wf_atoms(ispec)%get_psi(j, r)
          dens(i) = dens(i) + wf_atoms(ispec)%get_ref_charge(j)*psi(j)**2
        end do
      end if
      dens(i) = dens(i)*inv4pi
      if (.not. isder) cycle
      ddens(i, 1) = 0.0_dp
      dddens(i, 1, 1) = 0.0_dp
      if (r > tolerance) then
        do j = 1, wf_atoms(ispec)%get_nshells()
          tpsi = wf_atoms(ispec)%get_psi(j, r, order=1)
          ddens(i, 1) = ddens(i, 1) + wf_atoms(ispec)%get_ref_charge(j)*psi(j)*tpsi
         dddens(i, 1, 1) = dddens(i, 1, 1) + wf_atoms(ispec)%get_ref_charge(j) * &
           &               (tpsi**2 + psi(j)*wf_atoms(ispec)%get_psi(j, r, order=2))
        end do
      end if
      ddens(i, 1) = ddens(i, 1)*2.0_dp*inv4pi
      dddens(i, 1, 1) = dddens(i, 1, 1)*2.0_dp*inv4pi
    end do
    deallocate (psi)
  end subroutine dens1c

  pure subroutine dens2c(ispec1, ispec2, z1, z2, dens, ddens, dddens)
    integer, intent(in) :: ispec1, ispec2
    real(dp), intent(in) :: z1, z2
    real(dp), intent(out) :: dens(:)
    real(dp), intent(out), optional :: ddens(:, :), dddens(:, :, :)
    logical :: isder
    integer :: i, k, nrho, nz
    real(dp) :: rho, rho2, z12, z22, zmin, zmax, drho, dz, r1, r2, ir, ir2, tpsi
    real(dp), allocatable :: psi1(:), psi2(:)
    nrho = size(dens, 1)
    drho = min(wf_atoms(ispec1)%get_rcut(), wf_atoms(ispec2)%get_rcut())/real(nrho - 1, kind=dp)
    z12 = z1*z1
    z22 = z2*z2
    isder = present(ddens) .and. present(dddens)
    allocate (psi1(wf_atoms(ispec1)%get_nshells()), psi2(wf_atoms(ispec2)%get_nshells()))
    do i = 1, nrho
      rho = real(i - 1, kind=dp)*drho
      rho2 = rho*rho
      r1 = sqrt(rho2 + z12)
      r2 = sqrt(rho2 + z22)
      dens(i) = 0.0_dp
      if (r1 > tolerance) then
        do k = 1, wf_atoms(ispec1)%get_nshells()
          psi1(k) = wf_atoms(ispec1)%get_psi(k, r1)
          dens(i) = dens(i) + wf_atoms(ispec1)%get_ref_charge(k)*psi1(k)**2
        end do
      end if
      if (r2 > tolerance) then
        do k = 1, wf_atoms(ispec2)%get_nshells()
          psi2(k) = wf_atoms(ispec2)%get_psi(k, r2)
          dens(i) = dens(i) + wf_atoms(ispec2)%get_ref_charge(k)*psi2(k)**2
        end do
      end if
      dens(i) = dens(i)*inv4pi
      if (.not. isder) cycle
      ddens(i, :) = 0.0_dp
      dddens(i, :, :) = 0.0_dp
      if (r1 > tolerance) then
        ir = 1.0_dp/r1
        ir2 = ir*ir
        do k = 1, wf_atoms(ispec1)%get_nshells()
          tpsi = wf_atoms(ispec1)%get_psi(k, r1, order=1)
          ddens(i, 1) = ddens(i, 1) + wf_atoms(ispec1)%get_ref_charge(k)*psi1(k)
          ddens(i, 2) = ddens(i, 2) + wf_atoms(ispec1)%get_ref_charge(k)*psi1(k)
          dddens(i, 1, 1) = dddens(i, 1, 1) + wf_atoms(ispec1)%get_ref_charge(k)* &
            &               (tpsi**2*rho2 + psi1(k)*wf_atoms(ispec1)%get_psi(k, r1, order=2)*rho2 + psi1(k)*tpsi*z12*ir)
          dddens(i, 1, 2) = dddens(i, 1, 2) + wf_atoms(ispec1)%get_ref_charge(k)* &
            &               (tpsi**2 + psi1(k)*wf_atoms(ispec1)%get_psi(k, r1, order=2) - psi1(k)*tpsi*ir)
          dddens(i, 2, 2) = dddens(i, 2, 2) + wf_atoms(ispec1)%get_ref_charge(k)* &
            &               (tpsi**2*z12 + psi1(k)*wf_atoms(ispec1)%get_psi(k, r1, order=2)*z12 + psi1(k)*tpsi*rho2*ir)
        end do
        ddens(i, 1) = ddens(i, 1)*rho*ir
        ddens(i, 2) = ddens(i, 2)*z1*ir
        dddens(i, 1, 1) = dddens(i, 1, 1)*ir2
        dddens(i, 1, 2) = dddens(i, 1, 2)*ir2*rho*z1
        dddens(i, 2, 2) = dddens(i, 2, 2)*ir2
      end if
      if (r2 > tolerance) then
        ir = 1.0_dp/r2
        ir2 = ir*ir
        do k = 1, wf_atoms(ispec2)%get_nshells()
          tpsi = wf_atoms(ispec2)%get_psi(k, r2, order=1)
          ddens(i, 1) = ddens(i, 1) + wf_atoms(ispec2)%get_ref_charge(k)*psi2(k)
          ddens(i, 2) = ddens(i, 2) + wf_atoms(ispec2)%get_ref_charge(k)*psi2(k)
          dddens(i, 1, 1) = dddens(i, 1, 1) + wf_atoms(ispec2)%get_ref_charge(k)* &
            &               (tpsi**2*rho2 + psi2(k)*wf_atoms(ispec2)%get_psi(k, r2, order=2)*rho2 + psi2(k)*tpsi*z22*ir)
          dddens(i, 1, 2) = dddens(i, 1, 2) + wf_atoms(ispec2)%get_ref_charge(k)* &
            &               (tpsi**2 + psi2(k)*wf_atoms(ispec2)%get_psi(k, r2, order=2) - psi2(k)*tpsi*ir)
          dddens(i, 2, 2) = dddens(i, 2, 2) + wf_atoms(ispec2)%get_ref_charge(k)* &
            &               (tpsi**2*z22 + psi2(k)*wf_atoms(ispec2)%get_psi(k, r2, order=2)*z22 + psi2(k)*tpsi*rho2*ir)
        end do
        ddens(i, 1) = ddens(i, 1)*rho*ir
        ddens(i, 2) = ddens(i, 2)*z2*ir
        dddens(i, 1, 1) = dddens(i, 1, 1)*ir2
        dddens(i, 1, 2) = dddens(i, 1, 2)*ir2*rho*z2
        dddens(i, 2, 2) = dddens(i, 2, 2)*ir2
      end if
      ddens(i, :) = ddens(i, :)*2.0_dp*inv4pi
      dddens(i, 2, 1) = dddens(i, 1, 2)
      dddens(i, :, :) = dddens(i, :, :)*2.0_dp*inv4pi
    end do
    deallocate (psi1, psi2)
  end subroutine dens2c

end module wavefunctions
