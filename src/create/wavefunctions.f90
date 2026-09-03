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
    real(dp) :: rcut, rcut_max, q, qref, dr
    type(math_interp_t) :: fr
  contains
    private
    procedure :: end => wf_orbital_end
  end type wf_orbital_t

  type :: wf_atom_t
    private
    integer :: nz, nshells
    real(dp) :: rcut, dr
    type(wf_orbital_t), allocatable :: wfs(:)
  contains
    private
    procedure, public :: get_nz => wf_atom_nz
    procedure, public :: get_nshells => wf_atom_nshells
    procedure, public :: get_angular_momentum => wf_atom_angular_momentum
    procedure, public :: get_charge => wf_atom_charge
    procedure, public :: get_ref_charge => wf_atom_ref_charge
    procedure, public :: get_rcut => wf_atom_rcut
    procedure, public :: get_dr => wf_atom_dr
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
    real(dp) :: rcut, dr
    type(wf_orbital_t), allocatable :: wfs(:)
    allocate (wfs(nshells))
    rcut = 0.0_dp
    dr = 1000000.0_dp
    do i = 1, nshells
      wfs(i) = wf_new_orbital(fpaths(i))
      rcut = max(rcut, wfs(i)%rcut)
      dr = min(dr, wfs(i)%dr)
    end do
    wf_new_atom = wf_atom_t(nz=nz, nshells=nshells, rcut=rcut, dr=dr, wfs=wfs)
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
    call wf_normalize(fr)
    deallocate (r, psi)
    ! TODO: allow qref to be different
    wf_new_orbital = wf_orbital_t(l=l, rcut=rcut, rcut_max=rcut_max, q=q, qref=q, dr=dr, fr=fr)
  end function wf_new_orbital

  subroutine wf_normalize(fr)
    class(math_interp_t), intent(inout) :: fr
    integer :: i
    real(dp) :: normsq, r, r2, x, a, a2, b, b2, c, c2, d, d2, ab, ac, ad, bc, bd, cd
    normsq = 0.0_dp
    do i = 1, (fr%get_n() - 1)
      r = fr%get_x(i)
      x = fr%get_x(i + 1) - r
      a = fr%get_coef(1, i)
      b = fr%get_coef(2, i)
      c = fr%get_coef(3, i)
      d = fr%get_coef(4, i)
      r2 = r*r
      a2 = a*a
      b2 = b*b
      c2 = c*c
      d2 = d*d
      ab = a*b
      ac = a*c
      ad = a*d
      bc = b*c
      bd = b*d
      cd = c*d
      normsq = normsq + x * ( &
        &      a2*r2 + x * ( &
        &      ab*r2 + a2*r + x * ( &
        &      0.33333333333333333_dp*(b2*r2 + 2.0_dp*ac*r2 + 4.0_dp*ab*r + a2) + x * ( &
        &      0.5_dp*(ad*r2 + bc*r2 + b2*r + 2.0_dp*ac*r + ab) + x * ( &
        &      0.2_dp*(c2*r2 + 2.0_dp*bd*r2 + 4.0_dp*ad*r + 4.0_dp*bc*r + b2 + 2.0_dp*ac) + x * ( &
        &      0.33333333333333333_dp*(cd*r2 + c2*r + 2.0_dp*bd*r + ad + bc) + x * ( &
        &      0.14285714285714285_dp*(d2*r2 + 4.0_dp*cd*r + c2 + 2.0_dp*bd) + x * ( &
        &      0.25_dp*(d2*r + cd) + x * ( &
        &      0.11111111111111111_dp*d2)))))))))
    end do
    call fr%rescale(1.0_dp/sqrt(normsq))
  end subroutine wf_normalize

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

  pure real(dp) function wf_atom_dr(this, ish)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in), optional :: ish
    if (present(ish)) then
      wf_atom_dr = this%wfs(ish)%dr
    else
      wf_atom_dr = this%dr
    end if
  end function wf_atom_dr

  pure real(dp) function wf_atom_psi(this, ish, r, order)
    class(wf_atom_t), intent(in) :: this
    integer, intent(in) :: ish
    real(dp), intent(in) :: r
    integer, intent(in), optional :: order
    integer :: o
    !if ((r <= 0.0_dp) .or. (r >= this%wfs(ish)%rcut)) then
    !  wf_atom_psi = 0.0_dp
    !  return
    !end if
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

  pure subroutine dens1c(ispec, r, dens, ddens, dddens)
    integer, intent(in) :: ispec
    real(dp), intent(in) :: r
    real(dp), intent(out) :: dens
    real(dp), intent(out), optional :: ddens(:), dddens(:,:)
    integer :: issh, nssh
    real(dp) :: tpsi, tdpsi
    real(dp), allocatable :: psi(:)
    allocate (psi(wf_atoms(ispec)%get_nshells()))
    nssh = wf_atoms(ispec)%get_nshells()
    dens = 0.0_dp
    do issh = 1, nssh
      tpsi = wf_atoms(ispec)%get_psi(issh, r)
      psi(issh) = tpsi
      dens = dens + wf_atoms(ispec)%get_ref_charge(issh)*tpsi*tpsi
    end do
    dens = dens*inv4pi
    if (.not. present(ddens) .or. .not. present(dddens)) return
    ddens = 0.0_dp
    dddens = 0.0_dp
    if (r > tolerance) then
      do issh = 1, nssh
        tpsi = psi(issh)
        tdpsi = wf_atoms(ispec)%get_psi(issh, r, order=1)
        ddens = ddens + wf_atoms(ispec)%get_ref_charge(issh)*tpsi*tdpsi
        dddens = dddens + wf_atoms(ispec)%get_ref_charge(issh) * &
          &               (tdpsi*tdpsi + tpsi*wf_atoms(ispec)%get_psi(issh, r, order=2))
      end do
    end if
    ddens = ddens*2.0_dp*inv4pi
    dddens = dddens*2.0_dp*inv4pi
    deallocate (psi)
  end subroutine dens1c

  pure subroutine dens2c(ispec, jspec, rho, z1, z2, r1, r2, dens, ddens, dddens)
    integer, intent(in) :: ispec, jspec
    real(dp), intent(in) :: rho, z1, z2, r1, r2
    real(dp), intent(out) :: dens
    real(dp), intent(out), optional :: ddens(:), dddens(:, :)
    integer :: issh, jssh, nssh1, nssh2
    real(dp) :: ir, ir2, tpsi, tdpsi, rho2, z12, z22
    real(dp), allocatable :: psi1(:), psi2(:)
    nssh1 = wf_atoms(ispec)%get_nshells()
    nssh2 = wf_atoms(jspec)%get_nshells()
    rho2 = rho*rho
    z12 = z1*z1
    z22 = z2*z2
    allocate (psi1(nssh1), psi2(nssh2))
    dens = 0.0_dp
    do issh = 1, nssh1
      tpsi = wf_atoms(ispec)%get_psi(issh, r1)
      psi1(issh) = tpsi
      dens = dens + wf_atoms(ispec)%get_ref_charge(issh)*tpsi*tpsi
    end do
    do jssh = 1, nssh2
      tpsi = wf_atoms(jspec)%get_psi(jssh, r2)
      psi2(jssh) = tpsi
      dens = dens + wf_atoms(jspec)%get_ref_charge(jssh)*tpsi*tpsi
    end do
    dens = dens*inv4pi
    if (.not. present(ddens) .or. .not. present(dddens)) return
    ddens = 0.0_dp
    dddens = 0.0_dp
    if (r1 > tolerance) then
      ir = 1.0_dp/r1
      ir2 = ir*ir
      do issh = 1, nssh1
        tdpsi = wf_atoms(ispec)%get_psi(issh, r1, order=1)
        tpsi = psi1(issh)
        ddens(1) = ddens(1) + wf_atoms(ispec)%get_ref_charge(issh)*tpsi*tdpsi
        ddens(2) = ddens(2) + wf_atoms(ispec)%get_ref_charge(issh)*tpsi*tdpsi
        dddens(1, 1) = dddens(1, 1) + wf_atoms(ispec)%get_ref_charge(issh)* &
          &            (tdpsi*tdpsi*rho2 + tpsi*wf_atoms(ispec)%get_psi(issh, r1, order=2)*rho2 + tpsi*tdpsi*z12*ir)
        dddens(1, 2) = dddens(1, 2) + wf_atoms(ispec)%get_ref_charge(issh)* &
          &            (tdpsi*tdpsi + tpsi*wf_atoms(ispec)%get_psi(issh, r1, order=2) - tpsi*tdpsi*ir)
        dddens(2, 2) = dddens(2, 2) + wf_atoms(ispec)%get_ref_charge(issh)* &
          &            (tdpsi*tdpsi*z12 + tpsi*wf_atoms(ispec)%get_psi(issh, r1, order=2)*z12 + tpsi*tdpsi*rho2*ir)
      end do
      ddens(1) = ddens(1)*rho*ir
      ddens(2) = ddens(2)*z1*ir
      dddens(1, 1) = dddens(1, 1)*ir2
      dddens(1, 2) = dddens(1, 2)*ir2*rho*z1
      dddens(2, 2) = dddens(2, 2)*ir2
    end if
    if (r2 > tolerance) then
      ir = 1.0_dp/r2
      ir2 = ir*ir
      do jssh = 1, nssh2
        tdpsi = wf_atoms(jspec)%get_psi(jssh, r2, order=1)
        tpsi = psi2(jssh)
        ddens(1) = ddens(1) + wf_atoms(jspec)%get_ref_charge(jssh)*tpsi*tdpsi
        ddens(2) = ddens(2) + wf_atoms(jspec)%get_ref_charge(jssh)*tpsi*tdpsi
        dddens(1, 1) = dddens(1, 1) + wf_atoms(jspec)%get_ref_charge(jssh)* &
          &            (tdpsi*tdpsi*rho2 + tpsi*wf_atoms(jspec)%get_psi(jssh, r2, order=2)*rho2 + tpsi*tdpsi*z22*ir)
        dddens(1, 2) = dddens(1, 2) + wf_atoms(jspec)%get_ref_charge(jssh)* &
          &            (tdpsi*tdpsi + tpsi*wf_atoms(jspec)%get_psi(jssh, r2, order=2) - tpsi*tdpsi*ir)
        dddens(2, 2) = dddens(2, 2) + wf_atoms(jspec)%get_ref_charge(jssh)* &
          &            (tdpsi*tdpsi*z22 + tpsi*wf_atoms(jspec)%get_psi(jssh, r2, order=2)*z22 + tpsi*tdpsi*rho2*ir)
      end do
      ddens(1) = ddens(1)*rho*ir
      ddens(2) = ddens(2)*z2*ir
      dddens(1, 1) = dddens(1, 1)*ir2
      dddens(1, 2) = dddens(1, 2)*ir2*rho*z2
      dddens(2, 2) = dddens(2, 2)*ir2
    end if
    ddens = ddens*2.0_dp*inv4pi
    dddens = dddens*2.0_dp*inv4pi
    dddens(2, 1) = dddens(1, 2)
    deallocate (psi1, psi2)
  end subroutine dens2c

end module wavefunctions
