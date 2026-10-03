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
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! xc.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of Exc, Vxc
!       and their partial derivatives using libXC
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module xc
  use, intrinsic :: iso_fortran_env, only: int64, dp => real64, stdout => output_unit, stderr => error_unit
  use :: constants, only: abohr3, abohr4, abohr5, abohr6, abohr8, abohr13, hartree, tolerance
  use :: xc_f03_lib_m, only: xc_f03_version, xc_f03_func_init_flags, xc_f03_func_get_info, xc_f03_func_info_get_family, &
    &                        xc_f03_func_end, xc_f03_lda_exc_vxc, xc_f03_lda_fxc, xc_f03_lda_kxc, &
    &                        xc_f03_gga_exc_vxc, xc_f03_gga_fxc, xc_f03_gga_kxc, xc_f03_func_t, xc_f03_func_info_t, &
    &                        XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FLAGS_ON_HOST
  implicit none
  private
  public :: xc_isgga, xc_calc, xc_init, xc_end, xc_get_nxc, xc_get_iexc, xc_get_weight

  interface xc_calc
    procedure :: xc_calc_lda
    procedure :: xc_calc_gga
  end interface xc_calc

  logical :: xc_gga  ! We force either all LDA or all GGA
  integer :: nxc
  integer, allocatable :: xc_iexcs(:), xc_family(:)
  real(dp), allocatable :: xc_weights(:)
  type(xc_f03_func_t), allocatable :: xc_func(:)
  type(xc_f03_func_info_t), allocatable :: xc_info(:)

contains

  pure logical function xc_isgga()
    xc_isgga = xc_gga
  end function xc_isgga

  pure integer function xc_get_nxc()
    xc_get_nxc = nxc
  end function xc_get_nxc

  pure integer function xc_get_iexc(iexc)
    integer, intent(in) :: iexc
    xc_get_iexc = xc_iexcs(iexc)
  end function xc_get_iexc

  pure real(dp) function xc_get_weight(iexc)
    integer, intent(in) :: iexc
    xc_get_weight = xc_weights(iexc)
  end function xc_get_weight

  subroutine xc_init(iexcs, weights, verbose)
    integer, intent(in) :: iexcs(:)
    real(dp), intent(in), optional :: weights(:)
    logical, intent(in), optional :: verbose
    integer :: i, vmajor, vminor, vmicro
    logical, allocatable :: isgga(:)
    real(dp), allocatable :: ws(:)
    nxc = size(iexcs)
    if (allocated(xc_iexcs)) deallocate(xc_iexcs)
    if (allocated(xc_family)) deallocate(xc_family)
    if (allocated(xc_func)) deallocate(xc_func)
    if (allocated(xc_info)) deallocate(xc_info)
    allocate(xc_iexcs(nxc), xc_family(nxc), xc_weights(nxc), xc_func(nxc), xc_info(nxc), isgga(nxc), ws(nxc))

    if (present(weights)) then
      if (size(weights) /= nxc) then
        write (stderr, "(a)") "[ERROR]: size of weights is not consistent with the number of XCs"
        error stop
      end if
      ws = weights
    else
      ws = 1.0_dp
    end if

    do i = 1, nxc
      call xc_f03_func_init_flags(xc_func(i), iexcs(i), XC_UNPOLARIZED, XC_FLAGS_ON_HOST)
      xc_iexcs(i) = iexcs(i)
      xc_weights(i) = ws(i)
      xc_info(i) = xc_f03_func_get_info(xc_func(i))
      xc_family(i) = xc_f03_func_info_get_family(xc_info(i))
      select case (xc_family(i))
      case (XC_FAMILY_LDA)
        isgga(i) = .false.
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        isgga(i) = .true.
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
        error stop
      end select
    end do
    if (any(isgga) .and. .not. all(isgga)) then
      write (stderr, "(a)") "[ERROR]: all functionals must be either LDA or GGA"
      error stop
    end if
    xc_gga = any(isgga)
    if (present(verbose)) then
      if (verbose) then
        call xc_f03_version(vmajor, vminor, vmicro)
        write (stdout, "('  Using libXC v',i1,'.',i1,'.',i1)") vmajor, vminor, vmicro
      end if
    end if
  end subroutine xc_init

  subroutine xc_end()
    integer :: i
    do i = 1, nxc
      call xc_f03_func_end(xc_func(i))
    end do
    deallocate(xc_iexcs, xc_family, xc_weights, xc_func, xc_info)
  end subroutine xc_end

  subroutine xc_calc_lda(dens, exc, vxc, dexc, dvxc, d2exc, d2vxc)
    real(dp), intent(in) :: dens
    real(dp), intent(out) :: exc, vxc, dexc, dvxc
    real(dp), intent(out), optional :: d2exc, d2vxc
    integer :: i
    real(dp) :: rho(1), irho(1), sigma(1), e(1), vrho(1), vsigma(1), &
      &         v2rho2(1), v2rhosigma(1), v2sigma2(1), &
      &         v3rho3(1), v3rho2sigma(1), v3rhosigma2(1), v3sigma3(1)
    logical :: need2
    need2 = present(d2exc) .and. present(d2vxc)
    exc = 0.0_dp
    vxc = 0.0_dp
    dexc = 0.0_dp
    dvxc = 0.0_dp
    if (need2) d2vxc = 0.0_dp
    rho(1) = dens*abohr3
    irho(1) = 1.0_dp/max(tolerance, rho(1))
    do i = 1, nxc
      select case (xc_family(i))
      case (XC_FAMILY_LDA)
        call xc_f03_lda_exc_vxc(xc_func(i), 1_int64, rho, e, vrho)
        call xc_f03_lda_fxc(xc_func(i), 1_int64, rho, v2rho2)
        if (need2) call xc_f03_lda_kxc(xc_func(i), 1_int64, rho, v3rho3)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        sigma(1) = 0.0_dp
        call xc_f03_gga_exc_vxc(xc_func(i), 1_int64, rho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func(i), 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
        if (need2) call xc_f03_gga_kxc(xc_func(i), 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
        stop
      end select
      exc = exc + xc_weights(i)*e(1)
      vxc = vxc + xc_weights(i)*vrho(1)
      dexc = dexc + xc_weights(i)*irho(1)*(vrho(1) - e(1))
      dvxc = dvxc + xc_weights(i)*v2rho2(1)
      if (need2) d2vxc = d2vxc + xc_weights(i)*v3rho3(1)
    end do
    exc = exc*hartree
    vxc = vxc*hartree
    dexc = dexc*hartree*abohr3
    dvxc = dvxc*hartree*abohr3
    if (need2) then
      d2exc = irho(1)*(dvxc - 2.0_dp*dexc)*abohr3
      d2vxc = d2vxc*hartree*abohr6
    end if
  end subroutine xc_calc_lda

  subroutine xc_calc_gga(dens, lapl, grad, hess, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross)
    real(dp), intent(in) :: dens, lapl
    real(dp), intent(in) :: grad(:), hess(:, :)
    real(dp), intent(out) :: exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross
    integer :: ndim, i, j
    real(dp) :: rho(1), irho(1), sigma(1), laplacian(1), crossed(1), e(1), &
      &         vrho(1), vsigma(1), v2rho2(1), v2rhosigma(1), v2sigma2(1), &
      &         v3rho3(1), v3rho2sigma(1), v3rhosigma2(1), v3sigma3(1)
    exc = 0.0_dp
    vxc = 0.0_dp
    dexcrho = 0.0_dp
    dvxcrho = 0.0_dp
    dexcsigma = 0.0_dp
    dvxcsigma = 0.0_dp
    dvxclapl = 0.0_dp
    dvxccross = 0.0_dp
    ndim = size(grad)
    rho(1) = abohr3*dens
    irho(1) = 1.0_dp/max(tolerance, rho(1))
    laplacian(1) = abohr5*lapl
    sigma(1) = 0.0_dp
    crossed(1) = 0.0_dp
    do i = 1, ndim
      sigma(1) = sigma(1) + grad(i)*grad(i)
      do j = 1, ndim
        crossed(1) = crossed(1) + grad(i)*hess(j, i)*grad(j)
      end do
    end do
    sigma(1) = abohr8*sigma(1)
    crossed(1) = abohr13*crossed(1)
    do i = 1, nxc
      select case (xc_family(i))
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_f03_gga_exc_vxc(xc_func(i), 1_int64, rho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func(i), 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
        call xc_f03_gga_kxc(xc_func(i), 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
        exc = exc + xc_weights(i)*e(1)
        vxc = vxc + xc_weights(i)*(vrho(1) - &
          &   2.0_dp*(vsigma(1)*laplacian(1) + v2rhosigma(1)*sigma(1) + 2.0_dp*v2sigma2(1)*crossed(1)))
        dexcrho = dexcrho + xc_weights(i)*irho(1)*(vrho(1) - e(1))
        dexcsigma = dexcsigma + xc_weights(i)*irho(1)*vsigma(1)
        dvxcrho = dvxcrho + xc_weights(i)*(v2rho2(1) - &
          &       2.0_dp*(v2rhosigma(1)*laplacian(1) + v3rho2sigma(1)*sigma(1) + 2.0_dp*v3rhosigma2(1)*crossed(1)))
        dvxcsigma = dvxcsigma - xc_weights(i)*(v2rhosigma(1) - &
          &         2.0_dp*(v2sigma2(1)*laplacian(1) + v3rhosigma2(1)*sigma(1) + 2.0_dp*v3sigma3(1)*crossed(1)))
        dvxclapl = dvxclapl - xc_weights(i)*2.0_dp*vsigma(1)
        dvxccross = dvxccross - xc_weights(i)*4.0_dp*v2sigma2(1)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not GGA"
        error stop
      end select
    end do
    exc = exc*hartree
    vxc = vxc*hartree
    dexcrho = dexcrho*hartree*abohr3
    dexcsigma = dexcsigma*hartree*abohr8
    dvxcrho = dvxcrho*hartree*abohr3
    dvxcsigma = dvxcsigma*hartree*abohr8
    dvxclapl = dvxclapl*hartree*abohr5
    dvxccross = dvxccross*hartree*abohr13
  end subroutine xc_calc_gga
end module xc
