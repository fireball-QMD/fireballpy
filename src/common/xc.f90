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
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! xc.f90
! Program Description
! ==============================================================================
!       This module calculates handles the computation of Exc and Vxc
!       using libXC
! ==============================================================================
! Code written by:
! Cyra Roldán Piñero
! ==============================================================================
module xc
  use, intrinsic :: iso_fortran_env, only: int64, dp => real64, stdout => output_unit, stderr => error_unit
  use :: constants, only:abohr3, abohr4, abohr5, abohr8, abohr13, hartree, tolerance
  use :: xc_f03_lib_m, only:xc_f03_version, xc_f03_func_init_flags, xc_f03_func_get_info, xc_f03_func_info_get_family, &
    &                       xc_f03_func_end, xc_f03_lda_exc_vxc, xc_f03_lda_fxc, xc_f03_gga_exc_vxc, &
    &                       xc_f03_gga_fxc, xc_f03_gga_kxc, xc_f03_lda_kxc, xc_f03_func_t, xc_f03_func_info_t, &
    &                       XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FLAGS_ON_HOST
  implicit none
  private
  public :: xc_isgga, xc_calc, xc_init, xc_end

  interface xc_calc
    procedure :: xc_calc_lda
    procedure :: xc_calc_gga
  end interface xc_calc

  logical :: xc_sep, xc_gga
  integer :: xc_family1, xc_family2
  type(xc_f03_func_t) :: xc_func1, xc_func2
  type(xc_f03_func_info_t) :: xc_info1, xc_info2

contains

  subroutine xc_init(iexc1, iexc2, verbose)
    integer, intent(in) :: iexc1
    integer, intent(in), optional :: iexc2
    logical, intent(in), optional :: verbose
    integer :: vmajor, vminor, vmicro
    logical :: l1, l2
    xc_sep = present(iexc2)
    call xc_f03_func_init_flags(xc_func1, iexc1, XC_UNPOLARIZED, XC_FLAGS_ON_HOST)
    xc_info1 = xc_f03_func_get_info(xc_func1)
    xc_family1 = xc_f03_func_info_get_family(xc_info1)
    select case (xc_family1)
    case (XC_FAMILY_LDA)
      l1 = .false.
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      l1 = .true.
    case default
      write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
      stop
    end select
    l2 = .false.
    if (xc_sep) then
      call xc_f03_func_init_flags(xc_func2, iexc2, XC_UNPOLARIZED, XC_FLAGS_ON_HOST)
      xc_info2 = xc_f03_func_get_info(xc_func2)
      xc_family2 = xc_f03_func_info_get_family(xc_info2)
      select case (xc_family2)
      case (XC_FAMILY_LDA)
        l2 = .false.
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        l2 = .true.
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
        stop
      end select
    end if
    xc_gga = l1 .or. l2
    call xc_f03_version(vmajor, vminor, vmicro)
    if (present(verbose)) then
      if (verbose) write (stdout, "('  Using libXC v',i1,'.',i1,'.',i1)") vmajor, vminor, vmicro
    end if
  end subroutine xc_init

  subroutine xc_end()
    call xc_f03_func_end(xc_func1)
    if (xc_sep) call xc_f03_func_end(xc_func2)
  end subroutine xc_end

  pure logical function xc_isgga()
    xc_isgga = xc_gga
  end function xc_isgga

  ! d2exc y d2vxc son opcionales: d2exc = d2(exc)/d(dens)2 y d2vxc = d2(vxc)/d(dens)2.
  ! Se obtienen de la tercera derivada de la energia (kxc), asi que libXC debe
  ! estar compilado con MAXORDER >= 3.
  subroutine xc_calc_lda(dens, exc, vxc, dexc, dvxc, d2exc, d2vxc)
    real(dp), intent(in) :: dens
    real(dp), intent(out) :: exc, vxc, dexc, dvxc
    real(dp), intent(out), optional :: d2exc, d2vxc
    real(dp) :: rho(1), irho(1), sigma(1), e(1), vrho(1), v2rho2(1), v3rho3(1), &
      &         vsigma(1), v2rhosigma(1), v2sigma2(1), &
      &         v3rho2sigma(1), v3rhosigma2(1), v3sigma3(1)
    real(dp) :: ad2exc, ad2vxc
    logical :: need2
    need2 = present(d2exc) .or. present(d2vxc)
    exc = 0.0_dp
    vxc = 0.0_dp
    dexc = 0.0_dp
    dvxc = 0.0_dp
    ad2exc = 0.0_dp
    ad2vxc = 0.0_dp
    rho(1) = dens*abohr3
    irho(1) = 1.0_dp/max(tolerance, rho(1))
    select case (xc_family1)
    case (XC_FAMILY_LDA)
      call xc_f03_lda_exc_vxc(xc_func1, 1_int64, rho, e, vrho)
      call xc_f03_lda_fxc(xc_func1, 1_int64, rho, v2rho2)
      if (need2) call xc_f03_lda_kxc(xc_func1, 1_int64, rho, v3rho3)
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      sigma(1) = 0.0_dp
      call xc_f03_gga_exc_vxc(xc_func1, 1_int64, rho, sigma, e, vrho, vsigma)
      call xc_f03_gga_fxc(xc_func1, 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      if (need2) call xc_f03_gga_kxc(xc_func1, 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    case default
      write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
      stop
    end select
    exc = exc + e(1)
    vxc = vxc + vrho(1)
    dexc = dexc + irho(1)*(vrho(1) - e(1))
    dvxc = dvxc + v2rho2(1)
    ! v = e + n de/dn  =>  d2e/dn2 = (dv/dn - 2 de/dn)/n ;  d2v/dn2 = v3rho3
    if (need2) ad2exc = ad2exc + irho(1)*(v2rho2(1) - 2.0_dp*irho(1)*(vrho(1) - e(1)))
    if (need2) ad2vxc = ad2vxc + v3rho3(1)
    if (xc_sep) then
      select case (xc_family2)
      case (XC_FAMILY_LDA)
        call xc_f03_lda_exc_vxc(xc_func2, 1_int64, rho, e, vrho)
        call xc_f03_lda_fxc(xc_func2, 1_int64, rho, v2rho2)
        if (need2) call xc_f03_lda_kxc(xc_func2, 1_int64, rho, v3rho3)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        sigma(1) = 0.0_dp
        call xc_f03_gga_exc_vxc(xc_func2, 1_int64, rho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func2, 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
        if (need2) call xc_f03_gga_kxc(xc_func2, 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
        stop
      end select
      exc = exc + e(1)
      vxc = vxc + vrho(1)
      dexc = dexc + irho(1)*(vrho(1) - e(1))
      dvxc = dvxc + v2rho2(1)
      if (need2) ad2exc = ad2exc + irho(1)*(v2rho2(1) - 2.0_dp*irho(1)*(vrho(1) - e(1)))
      if (need2) ad2vxc = ad2vxc + v3rho3(1)
    end if
    exc = exc*hartree
    vxc = vxc*hartree
    dexc = dexc*hartree*abohr3
    dvxc = dvxc*hartree*abohr3
    if (present(d2exc)) d2exc = ad2exc*hartree*abohr3*abohr3
    if (present(d2vxc)) d2vxc = ad2vxc*hartree*abohr3*abohr3
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
    select case (xc_family1)
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call xc_f03_gga_exc_vxc(xc_func1, 1_int64, rho, sigma, e, vrho, vsigma)
      call xc_f03_gga_fxc(xc_func1, 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      call xc_f03_gga_kxc(xc_func1, 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      exc = exc + e(1)
      vxc = vxc + vrho(1) - 2.0_dp*(vsigma(1)*laplacian(1) + v2rhosigma(1)*sigma(1) + 2.0_dp*v2sigma2(1)*crossed(1))
      dexcrho = dexcrho + irho(1)*(vrho(1) - e(1))
      dexcsigma = dexcsigma + irho(1)*vsigma(1)
      dvxcrho = dvxcrho + v2rho2(1) - &
        &       2.0_dp*(v2rhosigma(1)*laplacian(1) + v3rho2sigma(1)*sigma(1) + 2.0_dp*v3rhosigma2(1)*crossed(1))
      dvxcsigma = dvxcsigma - v2rhosigma(1) - &
        &         2.0_dp*(v2sigma2(1)*laplacian(1) + v3rhosigma2(1)*sigma(1) + 2.0_dp*v3sigma3(1)*crossed(1))
      dvxclapl = dvxclapl - 2.0_dp*vsigma(1)
      dvxccross = dvxccross - 4.0_dp*v2sigma2(1)
    case default
      write (stderr, "(a)") "[ERROR]: selected functional is not GGA"
      stop
    end select

    if (xc_sep) then
      select case (xc_family2)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_f03_gga_exc_vxc(xc_func2, 1_int64, rho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func2, 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
        call xc_f03_gga_kxc(xc_func2, 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
        exc = exc + e(1)
        vxc = vxc + vrho(1) - 2.0_dp*(vsigma(1)*laplacian(1) + v2rhosigma(1)*sigma(1) + 2.0_dp*v2sigma2(1)*crossed(1))
        dexcrho = dexcrho + irho(1)*(vrho(1) - e(1))
        dexcsigma = dexcsigma + irho(1)*vsigma(1)
        dvxcrho = dvxcrho + v2rho2(1) - &
          &       2.0_dp*(v2rhosigma(1)*laplacian(1) + v3rho2sigma(1)*sigma(1) + 2.0_dp*v3rhosigma2(1)*crossed(1))
        dvxcsigma = dvxcsigma - v2rhosigma(1) - &
          &         2.0_dp*(v2sigma2(1)*laplacian(1) + v3rhosigma2(1)*sigma(1) + 2.0_dp*v3sigma3(1)*crossed(1))
        dvxclapl = dvxclapl - 2.0_dp*vsigma(1)
        dvxccross = dvxccross - 4.0_dp*v2sigma2(1)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not GGA"
        stop
      end select
    end if

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
