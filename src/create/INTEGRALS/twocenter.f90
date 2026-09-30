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
! Brigham Young University - Hao Wang Lawrence
! Livermore National Laboratory - Kurt Glaesemann
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
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! onecenter.f90
! Program Description
! ==============================================================================
!       This module calculates the two-center integrals.
! ==============================================================================
! Code written by:
! Cyra Roldan Pinero
! ==============================================================================
module twocenter
  use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit, stdout => output_unit
  use :: constants, only: sqinv4pi, twopi, invsq2, tolerance, fname_len, inv4pi
  use :: utils, only: utils_open, utils_itos
  use :: indices, only: indices_twocenter_set, INDICES_TWOCENTER, INDICES_TWOCENTER_DIPX, &
    &                   INDICES_TWOCENTER_DIPY, INDICES_TWOCENTER_COULOMB, INDICES_TWOCENTER_SPH, &
    &                   INDICES_NAMES_LEN
  use :: xc, only: xc_calc, xc_isgga
  use :: wavefunctions, only: wf_atoms, wf_dens
  use :: potentials, only: pot_atoms
  use :: pseudopotentials, only: pp_atoms
  implicit none
  private
  public :: twocenter_calc

  integer, parameter, public :: TWOCENTER_NPOINTS_D = 107
  integer, parameter, public :: TWOCENTER_NPOINTS_Z = 213
  integer, parameter, public :: TWOCENTER_NPOINTS_RHO = 107

  integer, parameter :: TWOCENTER_NUM_INTERACTIONS = 13
  integer, parameter, public :: TWOCENTER_DENS_ATOM      = ishft(1, 0)
  integer, parameter, public :: TWOCENTER_DENS_ONTOP     = ishft(1, 1)
  integer, parameter, public :: TWOCENTER_OVERLAP        = ishft(1, 2)
  integer, parameter, public :: TWOCENTER_VNEUTRAL_ATOM  = ishft(1, 3)
  integer, parameter, public :: TWOCENTER_VNN_ATOM       = ishft(1, 4)
  integer, parameter, public :: TWOCENTER_VNEUTRAL_ONTOP = ishft(1, 5)
  integer, parameter, public :: TWOCENTER_VNN_ONTOP      = ishft(1, 6)
  integer, parameter, public :: TWOCENTER_VPP            = ishft(1, 7)
  integer, parameter, public :: TWOCENTER_VXC            = ishft(1, 8)
  integer, parameter, public :: TWOCENTER_DIP_Z          = ishft(1, 9)
  integer, parameter, public :: TWOCENTER_DIP_X          = ishft(1, 10)
  integer, parameter, public :: TWOCENTER_DIP_Y          = ishft(1, 11)
  integer, parameter, public :: TWOCENTER_COULOMB        = ishft(1, 12)

contains

  pure logical function twocenter_is_long_range(interaction)
    integer, intent(in) :: interaction
    select case (interaction)
    case (TWOCENTER_COULOMB, TWOCENTER_VNN_ATOM)
      twocenter_is_long_range = .true.
    case default
      twocenter_is_long_range = .false.
    end select
  end function twocenter_is_long_range

  pure integer function twocenter_get_nints(interaction, ispec, jspec)
    integer, intent(in) :: interaction, ispec, jspec
    select case (interaction)
    case (TWOCENTER_DENS_ATOM)
      twocenter_get_nints = wf_atoms(jspec)%get_nshells()
    case (TWOCENTER_VNN_ATOM)
      twocenter_get_nints = pot_atoms(jspec)%get_nshells()
    case (TWOCENTER_DENS_ONTOP)
      twocenter_get_nints = wf_atoms(ispec)%get_nshells() + wf_atoms(jspec)%get_nshells()
    case (TWOCENTER_VNN_ONTOP)
      twocenter_get_nints = pot_atoms(ispec)%get_nshells() + pot_atoms(jspec)%get_nshells()
    case (TWOCENTER_VXC)
      twocenter_get_nints = wf_atoms(ispec)%get_nshells() + wf_atoms(jspec)%get_nshells() + 1
    case (TWOCENTER_VNEUTRAL_ONTOP)
      twocenter_get_nints = 2
    case default
      twocenter_get_nints = 1
    end select
  end function twocenter_get_nints

  subroutine twocenter_get_fnames(interaction, ispec, jspec, issph, fnames)
    integer, intent(in) :: interaction, ispec, jspec
    logical, intent(in) :: issph
    character(fname_len), allocatable, intent(out) :: fnames(:)
    integer :: i, nints
    character(:), allocatable :: csph, tail
    character(:), allocatable :: roots(:)
    if (issph) then
      csph = "S"
    else
      csph = ""
    end if
    select case (interaction)
    case(TWOCENTER_DENS_ATOM)
      roots = ["den"//csph//"_atom"]
    case(TWOCENTER_DENS_ONTOP)
      roots = ["den"//csph//"_ontopl", "den"//csph//"_ontopr"]
    case(TWOCENTER_OVERLAP)
      roots = ["overlap"//csph]
    case(TWOCENTER_VNEUTRAL_ATOM)
      roots = ["vna"//csph//"_atom_00"]
    case(TWOCENTER_VNN_ATOM)
      roots = ["vna"//csph//"_atom"]
    case(TWOCENTER_VNEUTRAL_ONTOP)
      roots = ["vna"//csph//"_ontopl_00", "vna"//csph//"_ontopr_00"]
    case(TWOCENTER_VNN_ONTOP)
      roots = ["vna"//csph//"_ontopl", "vna"//csph//"_ontopr"]
    case(TWOCENTER_VPP)
      roots = ["vnl"//csph]
    case(TWOCENTER_VXC)
      roots = ["vxc"//csph//"_refchs", "vxc"//csph//"_ontopl", "vxc"//csph//"_ontopr"]
    case(TWOCENTER_DIP_Z)
      roots = ["dipole"//csph//"_z"]
    case(TWOCENTER_DIP_X)
      roots = ["dipole"//csph//"_x"]
    case(TWOCENTER_DIP_Y)
      roots = ["dipole"//csph//"_y"]
    case(TWOCENTER_COULOMB)
      roots = ["coulomb"//csph]
    case default
      error stop
    end select
    tail = "."//utils_itos(wf_atoms(ispec)%get_nz(), lpad=2)//"."//utils_itos(wf_atoms(jspec)%get_nz(), lpad=2)//".dat"

    nints = twocenter_get_nints(interaction, ispec, jspec)
    if (allocated(fnames)) deallocate(fnames)
    allocate(fnames(nints))
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNN_ATOM)
      fnames = [(roots(1)//"_"//utils_itos(i, lpad=2)//tail, i = 1, wf_atoms(jspec)%get_nshells())]
    case (TWOCENTER_DENS_ONTOP, TWOCENTER_VNN_ONTOP)
      fnames(:wf_atoms(ispec)%get_nshells()) = [(roots(1)//"_"//utils_itos(i, lpad=2)//tail, &
        &                                        i = 1, wf_atoms(ispec)%get_nshells())]
      fnames(wf_atoms(ispec)%get_nshells()+1:) = [(roots(2)//"_"//utils_itos(i, lpad=2)//tail, &
        &                                          i = 1, wf_atoms(jspec)%get_nshells())]
    case (TWOCENTER_VXC)
      fnames(1) = roots(1)//tail
      fnames(2:wf_atoms(ispec)%get_nshells()+1) = [(roots(2)//"_"//utils_itos(i, lpad=2)//tail, &
        &                                           i = 1, wf_atoms(ispec)%get_nshells())]
      fnames(wf_atoms(ispec)%get_nshells()+2:) = [(roots(3)//"_"//utils_itos(i, lpad=2)//tail, &
        &                                          i = 1, wf_atoms(jspec)%get_nshells())]
    case default
      fnames = [(roots(i)//tail, i = 1, nints)]
    end select
  end subroutine twocenter_get_fnames

  subroutine twocenter_get_rcuts(interaction, ispec, jspec, rcut1, rcut2)
    integer, intent(in) :: interaction, ispec, jspec
    real(dp), intent(out) :: rcut1, rcut2
    rcut1 = wf_atoms(ispec)%get_rcut()
    select case (interaction)
    case (TWOCENTER_VPP)
      rcut2 = pp_atoms(jspec)%get_rcut()
    case default
      rcut2 = wf_atoms(jspec)%get_rcut()
    end select
  end subroutine twocenter_get_rcuts

  subroutine twocenter_get_names(interaction, ispec, jspec, issph, index_max, names)
    integer, intent(in) :: interaction, ispec, jspec
    logical, intent(in) :: issph
    integer, intent(out) :: index_max
    character(INDICES_NAMES_LEN), allocatable, intent(out) :: names(:)
    integer, allocatable :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:)
    call twocenter_get_indices(interaction, ispec, jspec, issph, index_max, s1, s2, l1, l2, m1, m2, names=names)
  end subroutine twocenter_get_names

  subroutine twocenter_get_indices(interaction, ispec, jspec, issph, index_max, s1, s2, l1, l2, m1, m2, names)
    integer, intent(in) :: interaction, ispec, jspec
    logical, intent(in) :: issph
    integer, intent(out) :: index_max
    integer, allocatable, intent(out) :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:)
    character(INDICES_NAMES_LEN), allocatable, intent(out), optional :: names(:)
    integer :: i
    integer, allocatable :: ls1(:), ls2(:)
    ls1 = [(wf_atoms(ispec)%get_angular_momentum(i), i = 1, wf_atoms(ispec)%get_nshells())]
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNEUTRAL_ATOM, TWOCENTER_VNN_ATOM)
      ls2 = [(wf_atoms(ispec)%get_angular_momentum(i), i = 1, wf_atoms(ispec)%get_nshells())]
    case (TWOCENTER_VPP)
      ls2 = [(pp_atoms(jspec)%get_angular_momentum(i), i = 1, pp_atoms(jspec)%get_nshells())]
    case default
      ls2 = [(wf_atoms(jspec)%get_angular_momentum(i), i = 1, wf_atoms(jspec)%get_nshells())]
    end select

    if (issph) then
      call indices_twocenter_set(INDICES_TWOCENTER_SPH, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
      return
    end if
    select case (interaction)
    case (TWOCENTER_DIP_X)
      call indices_twocenter_set(INDICES_TWOCENTER_DIPX, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    case (TWOCENTER_DIP_Y)
      call indices_twocenter_set(INDICES_TWOCENTER_DIPY, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    case (TWOCENTER_COULOMB)
      call indices_twocenter_set(INDICES_TWOCENTER_COULOMB, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    case default
      call indices_twocenter_set(INDICES_TWOCENTER, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    end select
  end subroutine twocenter_get_indices

  subroutine twocenter_calc(interactions, ispec, jspec, is_spherical)
    integer, intent(in) :: interactions, ispec, jspec
    logical, intent(in), optional :: is_spherical
    integer :: interaction, i, j
    logical :: issph, aux, exists
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    character(fname_len), allocatable :: fnames(:)
    real(dp), allocatable :: answer(:,:,:)

    issph = .false.
    if (present(is_spherical)) issph = is_spherical
    call twocenter_get_interactions(interactions, twocenter_interactions)
    do j = 1, TWOCENTER_NUM_INTERACTIONS
      if (.not. twocenter_interactions(j)) cycle
      interaction = ishft(1, j - 1)
      call twocenter_get_fnames(interaction, ispec, jspec, issph, fnames)
      exists = .true.
      do i = 1, size(fnames)
        inquire (file=fnames(i), exist=aux)
        exists = exists .and. aux
      end do
      if (exists) cycle
      do i = 1, size(fnames)
        write (stdout, "(2x,a)") "Computing "//trim(fnames(i))//"..."
      end do
      call twocenter_integrate(interaction, ispec, jspec, issph, answer)
      call twocenter_write(interaction, ispec, jspec, issph, fnames, answer)
      write (stdout, "(a)", advance="no") repeat(achar(8)//achar(13), size(fnames))
      do i = 1, size(fnames)
        write (stdout, "(2x,a)") "Computing "//trim(fnames(i))//"... Done!"
      end do
    end do
  end subroutine twocenter_calc

  subroutine twocenter_point_eval(interaction, ispec, jspec, d, rho, z1, z2, r1, r2, fofr)
    integer, intent(in) :: interaction, ispec, jspec
    real(dp), intent(in) :: d, rho, z1, z2, r1, r2
    real(dp), intent(out) :: fofr(:)
    integer :: i
    real(dp) :: tmp, dtmp, ddtmp, ir1, ir2, gdggdgr, gdggrgd, &
      &         dens, lapl, exc, vxc, dexcrho, dvxcrho, dexcsigma, dvxcsigma, dvxclapl, dvxccross, grader
    real(dp), allocatable :: grad(:), hess(:,:)

    select case (interaction)
    case (TWOCENTER_DENS_ATOM)
      fofr = inv4pi*[(wf_atoms(jspec)%get_psi(i, r2)**2, i = 1, wf_atoms(jspec)%get_nshells())]
    case (TWOCENTER_DENS_ONTOP)
      fofr(:wf_atoms(ispec)%get_nshells()) = inv4pi*[(wf_atoms(ispec)%get_psi(i, r1)**2, i = 1, wf_atoms(ispec)%get_nshells())]
      fofr(wf_atoms(ispec)%get_nshells()+1:) = inv4pi*[(wf_atoms(jspec)%get_psi(i, r2)**2, i = 1, wf_atoms(jspec)%get_nshells())]
    case (TWOCENTER_VNEUTRAL_ATOM)
      fofr(1) = pot_atoms(jspec)%get_vneutral(r2)
    case (TWOCENTER_VNN_ATOM)
      fofr = [(pot_atoms(jspec)%get_vnn(i, r2), i = 1, pot_atoms(jspec)%get_nshells())]
    case (TWOCENTER_VNEUTRAL_ONTOP)
      fofr = [pot_atoms(ispec)%get_vneutral(r1), pot_atoms(jspec)%get_vneutral(r2)]
    case (TWOCENTER_VNN_ONTOP)
      fofr(:pot_atoms(ispec)%get_nshells()) = [(pot_atoms(ispec)%get_vnn(i, r1), i = 1, pot_atoms(ispec)%get_nshells())]
      fofr(pot_atoms(ispec)%get_nshells()+1:) = [(pot_atoms(jspec)%get_vnn(i, r2), i = 1, pot_atoms(jspec)%get_nshells())]
    case (TWOCENTER_VXC)
      if (xc_isgga()) then
        call wf_dens(ispec, jspec, rho, z1, z2, r1, r2, dens, lapl, grad, hess)
        call xc_calc(dens, lapl, grad, hess, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross)
        if (r1 > tolerance) ir1 = 1.0_dp/r1
        if (r2 > tolerance) ir2 = 1.0_dp/r2
      else
        call wf_dens(ispec, jspec, rho, z1, z2, r1, r2, dens)
        call xc_calc(dens, exc, vxc, dexcrho, dvxcrho)
      end if
      fofr(1) = vxc
      do i = 1, wf_atoms(ispec)%get_nshells()
        tmp = sqinv4pi*wf_atoms(ispec)%get_psi(i, r1)
        fofr(i + 1) = dvxcrho*tmp*tmp
        if (xc_isgga() .and. r1 > tolerance) then
          dtmp = sqinv4pi*wf_atoms(ispec)%get_psi(i, r1, order=1)
          ddtmp = sqinv4pi*wf_atoms(ispec)%get_psi(i, r1, order=2)
          grader = (grad(1)*rho + grad(2)*z1)*ir1
          gdggdgr = ir1*(grad(1)*(hess(1, 1)*rho + hess(1, 2)*z1) + grad(2)*(hess(2, 1)*rho + hess(2, 2)*z1))
          gdggrgd = ir1*(grad(1)*grad(1) + grad(2)*grad(2) - grader*grader)
          fofr(i + 1) = fofr(i + 1) + &
            &  dvxcsigma*4.0_dp*tmp*dtmp*grader + &
            &  dvxclapl*2.0_dp*(dtmp*dtmp + tmp*ddtmp + 2.0_dp*ir1*tmp*dtmp) + &
            &  dvxccross*2.0_dp*(grader*grader*(dtmp*dtmp + tmp*ddtmp) + tmp*dtmp*(gdggrgd + 2.0_dp*gdggdgr))
        end if
      end do
      do i = 1, wf_atoms(jspec)%get_nshells()
        tmp = sqinv4pi*wf_atoms(jspec)%get_psi(i, r2)
        fofr(i + 1 + wf_atoms(ispec)%get_nshells()) = dvxcrho*tmp*tmp
        if (xc_isgga() .and. r2 > tolerance) then
          dtmp = sqinv4pi*wf_atoms(jspec)%get_psi(i, r2, order=1)
          ddtmp = sqinv4pi*wf_atoms(jspec)%get_psi(i, r2, order=2)
          grader = (grad(1)*rho + grad(2)*z2)*ir2
          gdggdgr = ir2*(grad(1)*(hess(1, 1)*rho + hess(1, 2)*z2) + grad(2)*(hess(2, 1)*rho + hess(2, 2)*z2))
          gdggrgd = ir2*(grad(1)*grad(1) + grad(2)*grad(2) - grader*grader)
          fofr(i + wf_atoms(ispec)%get_nshells() + 1) = fofr(i + wf_atoms(ispec)%get_nshells() + 1) + &
            &  dvxcsigma*4.0_dp*tmp*dtmp*grader + &
            &  dvxclapl*2.0_dp*(dtmp*dtmp + tmp*ddtmp + 2.0_dp*ir2*tmp*dtmp) + &
            &  dvxccross*2.0_dp*(grader*grader*(dtmp*dtmp + tmp*ddtmp) + tmp*dtmp*(gdggrgd + 2.0_dp*gdggdgr))
        end if
      end do
    case (TWOCENTER_DIP_Z)
      fofr(1) = z1 - 0.5_dp*d
    case (TWOCENTER_DIP_X, TWOCENTER_DIP_Y)
      fofr(1) = rho
    case default
      fofr(1) = 1.0_dp
    end select
  end subroutine twocenter_point_eval

  real(dp) function twocenter_get_psimult(interaction, ispec, jspec, s1, s2, l1, l2, m1, m2, rho, z1, z2, r1, r2, issph)
    integer, intent(in) :: interaction, ispec, jspec, s1, s2, l1, l2, m1, m2
    real(dp), intent(in) :: rho, z1, z2, r1, r2
    logical, intent(in) :: issph
    real(dp) :: psi1, psi2, cyl1, cyl2

    psi1 = wf_atoms(ispec)%get_psi(s1, r1)
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNEUTRAL_ATOM, TWOCENTER_VNN_ATOM)
      psi2 = wf_atoms(ispec)%get_psi(s2, r1)
    case (TWOCENTER_VPP)
      psi2 = pp_atoms(jspec)%get_vpp(s2, r2)
    case (TWOCENTER_COULOMB)
      psi2 = psi1*pot_atoms(jspec)%get_vnn(s2, r2)
    case default
      psi2 = wf_atoms(jspec)%get_psi(s2, r2)
    end select

    cyl1 = cylindrical_harmonic(interaction, l1, m1, rho, z1, r1)
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNEUTRAL_ATOM, TWOCENTER_VNN_ATOM)
      cyl2 = cylindrical_harmonic(interaction, l2, m2, rho, z1, r1)
    case default
      cyl2 = cylindrical_harmonic(interaction, l2, m2, rho, z2, r2)
    end select

    twocenter_get_psimult = psi1*psi2
    if (issph .and. (twocenter_get_psimult < 0)) twocenter_get_psimult = -twocenter_get_psimult
    twocenter_get_psimult = twocenter_get_psimult*cyl1*cyl2
  end function twocenter_get_psimult

  subroutine twocenter_integrate(interaction, ispec, jspec, issph, answer)
    integer, intent(in) :: interaction, ispec, jspec
    logical, intent(in) :: issph
    real(dp), allocatable, intent(out) :: answer(:,:,:)
    integer :: nints, igrid, iz, irho, iint, idx, nz, nrho, index_max
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    real(dp) :: d, dd, dmax, rcut1, rcut2, zmin, zmax, dz, drho, z1, z2, rho, rhomult, rhomax, zmult, factor, r1, r2, psimult
    integer, allocatable :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:)
    real(dp), allocatable :: fofr(:)

    nints = twocenter_get_nints(interaction, ispec, jspec)
    call twocenter_get_indices(interaction, ispec, jspec, issph, index_max, s1, s2, l1, l2, m1, m2)
    if (allocated(answer)) deallocate(answer)
    allocate(fofr(nints), answer(nints, index_max, TWOCENTER_NPOINTS_D))
    if (index_max == 0) return
    answer = 0.0_dp

    ! Prepare integral
    call twocenter_get_rcuts(interaction, ispec, jspec, rcut1, rcut2)
    dmax = rcut1 + rcut2
    dd = dmax/real(TWOCENTER_NPOINTS_D - 1, kind=dp)
    ! Step size should be the same for all cases, we later adjust the number of points
    dz = (rcut1 + rcut2)/real(TWOCENTER_NPOINTS_Z - 1, kind=dp)
    drho = max(rcut1, rcut2)/real(TWOCENTER_NPOINTS_RHO - 1, kind=dp)

    do igrid = 1, TWOCENTER_NPOINTS_D
      d = real(igrid - 1, kind=dp)*dd
      if (twocenter_is_long_range(interaction)) then
        zmin = -rcut1
        zmax = rcut1
        rhomax = rcut1
      else
        zmin = max(-rcut1, d - rcut2)
        zmax = min(rcut1, d + rcut2)
        rhomax = min(rcut1, rcut2)
      end if

      ! Get number of integration points and ensure is odd
      nz = int((zmax - zmin)/dz) + 1
      nrho = int(rhomax/drho) + 1
      if (iand(nz, 1) == 0) nz = nz + 1
      if (iand(nrho, 1) == 0) nrho = nrho + 1

      do iz = 1, nz
        z1 = zmin + real(iz - 1, kind=dp)*dz
        z2 = z1 - d
        zmult = 0.66666666666666666667_dp*dz
        if (iz == 1 .or. iz == nz) then
          zmult = 0.5_dp*zmult
        else if (iand(iz, 1) == 0) then
          zmult = 2.0_dp*zmult
        end if

        do irho = 2, nrho  ! We can ignore 0 because of factor
          rho = real(irho - 1, kind=dp)*drho
          rhomult = 0.66666666666666666667_dp*drho
          if (irho == nrho) then
            rhomult = 0.5_dp*rhomult
          else if (iand(irho, 1) == 0) then
            rhomult = 2.0_dp*rhomult
          end if
          factor = zmult*rhomult*twopi

          r1 = sqrt(z1*z1 + rho*rho)
          if (r1 >= wf_atoms(ispec)%get_rcut()) cycle  ! Placing the same condition for r2 is much more tricky
          r2 = sqrt(z2*z2 + rho*rho)
          call twocenter_point_eval(interaction, ispec, jspec, d, rho, z1, z2, r1, r2, fofr)
          if (all(fofr == 0.0_dp)) continue
          do idx = 1, index_max
            psimult = twocenter_get_psimult(interaction, ispec, jspec, s1(idx), s2(idx), l1(idx), l2(idx), &
              &                             m1(idx), m2(idx), rho, z1, z2, r1, r2, issph)
            answer(:, idx, igrid) = answer(:, idx, igrid) + fofr*psimult*factor*rho
          end do ! idx
        end do ! irho
      end do ! iz
    end do ! igrid
  end subroutine twocenter_integrate

  subroutine twocenter_write(interaction, ispec, jspec, issph, fnames, answer)
    integer, intent(in) :: interaction, ispec, jspec
    logical, intent(in) :: issph
    character(fname_len), intent(in) :: fnames(:)
    real(dp), intent(in) :: answer(:,:,:)
    integer :: i, ish, igrid, index, io, index_max
    real(dp) :: rcut1, rcut2
    character(INDICES_NAMES_LEN), allocatable :: names(:)

    call twocenter_get_rcuts(interaction, ispec, jspec, rcut1, rcut2)
    call twocenter_get_names(interaction, ispec, jspec, issph, index_max, names=names)
    do i = 1, size(fnames)
      io = utils_open(fnames(i), "w")
      write (io, "(14x,i2,2x,i2,40x,'! Atomic numbers')") wf_atoms(ispec)%get_nz(), wf_atoms(jspec)%get_nz()
      write (io, "(2x,f14.6,2x,f14.6,28x,'! Cutoff radii')") rcut1, rcut2
      write (io, "(2x,f14.6,2x,f14.6,2x,i8,18x,'! zmin, zmax, nz')") 0.0_dp, rcut1 + rcut2, TWOCENTER_NPOINTS_D
      if (interaction == TWOCENTER_VPP) then
        write (io, "(8x,i8,36x,'! npp')") pp_atoms(jspec)%get_nshells()
        write (io, "(1000ES16.8)") (pp_atoms(jspec)%get_cl(ish), ish = 1, pp_atoms(jspec)%get_nshells())
      end if
      if (interaction == TWOCENTER_VXC) then
        write (io, "(8x,i8,2x,i8,26x,'! Number of shells')") wf_atoms(ispec)%get_nshells(), wf_atoms(jspec)%get_nshells()
        write (io, "(1000ES16.8)", advance="no") (wf_atoms(ispec)%get_ref_charge(ish), ish = 1, wf_atoms(ispec)%get_nshells())
        write (io, "(a)") " ! Reference charges atom 1"
        write (io, "(1000ES16.8)", advance="no") (wf_atoms(jspec)%get_ref_charge(ish), ish = 1, wf_atoms(jspec)%get_nshells())
        write (io, "(a)") " ! Reference charges atom 2"
      end if
      write (io, "('!')", advance="no")
      if (index_max == 0) then
        write (io, "(a)") ""
        close (io)
        return
      end if
      do index = 1, index_max
        write (io, "(4x,a)", advance="no") trim(names(index))
      end do
      if (size(answer, 2) /= 0) then
        write (io, "(a)") ""
        do igrid = 1, TWOCENTER_NPOINTS_D
          write (io, "(1000ES16.8)") answer(i, :, igrid)
        end do
      end if
      close (io)
    end do
  end subroutine twocenter_write

  pure subroutine twocenter_get_interactions(interactions, twocenter_interactions)
    integer, intent(in) :: interactions
    logical, intent(out) :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    integer :: i, p
    twocenter_interactions = .false.
    p = ishft(1, TWOCENTER_NUM_INTERACTIONS - 1)
    do i = TWOCENTER_NUM_INTERACTIONS, 1, -1
      if (iand(interactions, p) == p) then
        twocenter_interactions(i) = .true.
      end if
      p = ishft(p, -1)
    end do
  end subroutine twocenter_get_interactions

  real(dp) function cylindrical_harmonic(interaction, l, m, rho, z, r)
  ! s(sigma)  = 1
  ! p_sigma   = z/r
  ! p_pi      = rho/r
  ! d_sigma   = (2*z**2-rho**2)/r**2
  ! d_pi      = rho*z/r**2
  ! d_delta   = rho**2/r**2
  ! f_sigma   = z*(2*z**2-3*rho**2)/r**3
  ! f_pi      = rho*(4*z**2-rho**2)/r**3
  ! f_delta   = z*rho**2/r**3
  ! f_phi     = rho**3/r**3
    integer, intent(in) :: interaction, l, m
    real(dp), intent(in) :: rho, z, r

    if (interaction == TWOCENTER_COULOMB) then
      cylindrical_harmonic = sqinv4pi
      return
    end if

    if (r < tolerance) then
      cylindrical_harmonic = 0.0_dp
      return
    end if

    select case (l)
    case (0)
      cylindrical_harmonic = 1.0_dp
    case (1)
      select case (m)
      case (-1, 1)
        cylindrical_harmonic = 1.224744871391589_dp*rho/r
      case (0)
        cylindrical_harmonic = 1.7320508075688772_dp*z/r
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (2)
      select case (m)
      case (-2, 2)
        cylindrical_harmonic = 1.3693063937629153_dp*(rho*rho)/(r*r)
      case (-1, 1)
        cylindrical_harmonic = 2.7386127875258306_dp*(rho*z)/(r*r)
      case (0)
        cylindrical_harmonic = 1.118033988749895_dp*(2.0_dp*z*z - rho*rho)/(r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (3)
      select case (m)
      case (-3, 3)
        cylindrical_harmonic = 1.479019945774904_dp*(rho*rho*rho)/(r*r*r)
      case (-2, 2)
        cylindrical_harmonic = 3.6228441865473595_dp*(rho*rho*z)/(r*r*r)
      case (-1, 1)
        cylindrical_harmonic = 1.14564392373896_dp*rho*(4.0_dp*z*z - rho*rho)/(r*r*r)
      case (0)
        cylindrical_harmonic = 1.3228756555322954_dp*z*(2.0_dp*z*z - 3.0_dp*rho*rho*rho)/(r*r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
      case default
        cylindrical_harmonic = 0.0_dp
    end select
    cylindrical_harmonic = cylindrical_harmonic * sqinv4pi

    ! Some signs and factors need to be adjusted for dipX and dipY
    if (interaction == TWOCENTER_DIP_X .or. interaction == TWOCENTER_DIP_Y) then
      if (iand(m, 1) == 0) cylindrical_harmonic = cylindrical_harmonic*invsq2
      if (interaction == TWOCENTER_DIP_Y) then
        if (l == 2 .and. m == 2 .or. &
          & l == 3 .and. (m == -3 .or. m == 2 .or. m == 3)) cylindrical_harmonic = -cylindrical_harmonic
      end if
    end if
  end function cylindrical_harmonic
end module twocenter
