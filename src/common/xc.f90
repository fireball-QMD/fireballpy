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
  use :: xc_f03_lib_m, only:xc_f03_version, xc_f03_func_init, xc_f03_func_get_info, xc_f03_func_info_get_family, &
    &                       xc_f03_func_end, xc_f03_lda_exc_vxc, xc_f03_lda_fxc, xc_f03_gga_exc_vxc, &
    &                       xc_f03_gga_fxc, xc_f03_gga_kxc, xc_f03_func_t, xc_f03_func_info_t, &
    &                       XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA
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

  subroutine xc_init(iexc1, iexc2)
    integer, intent(in) :: iexc1
    integer, intent(in), optional :: iexc2
    integer :: vmajor, vminor, vmicro
    logical :: l1, l2
    xc_sep = present(iexc2)
    call xc_f03_func_init(xc_func1, iexc1, XC_UNPOLARIZED)
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
      call xc_f03_func_init(xc_func2, iexc2, XC_UNPOLARIZED)
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
    write (stdout, "('  Using libXC v',i1,'.',i1,'.',i1)") vmajor, vminor, vmicro
  end subroutine xc_init

  subroutine xc_end()
    call xc_f03_func_end(xc_func1)
    if (xc_sep) call xc_f03_func_end(xc_func2)
  end subroutine xc_end

  pure logical function xc_isgga()
    xc_isgga = xc_gga
  end function xc_isgga

  subroutine xc_calc_lda(rho, exc, vxc, dexc, dvxc)
    real(dp), intent(in) :: rho(:)
    real(dp), intent(out) :: exc(:), vxc(:), dexc(:), dvxc(:)
    integer(int64) :: np
    real(dp), allocatable :: nrho(:), irho(:), sigma(:), e(:), vrho(:), v2rho2(:), &
      &                      vsigma(:), v2rhosigma(:), v2sigma2(:)
    np = size(rho, kind=int64)
    allocate (nrho(np), irho(np), e(np), vrho(np), v2rho2(np))
    nrho = rho*abohr3
    irho = 1.0_dp/max(tolerance, nrho)
    exc = 0.0_dp
    vxc = 0.0_dp
    dexc = 0.0_dp
    dvxc = 0.0_dp
    select case (xc_family1)
    case (XC_FAMILY_LDA)
      call xc_f03_lda_exc_vxc(xc_func1, np, nrho, e, vrho)
      call xc_f03_lda_fxc(xc_func1, np, nrho, v2rho2)
      exc = exc + e
      vxc = vxc + vrho
      dexc = dexc + irho*(vrho - e)
      dvxc = dvxc + v2rho2
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      allocate (sigma(np), vsigma(np), v2rhosigma(np), v2sigma2(np))
      sigma = 0.0_dp
      call xc_f03_gga_exc_vxc(xc_func1, np, nrho, sigma, e, vrho, vsigma)
      call xc_f03_gga_fxc(xc_func1, np, nrho, sigma, v2rho2, v2rhosigma, v2sigma2)
      exc = exc + e
      vxc = vxc + vrho
      dexc = dexc + irho*(vrho - e)
      dvxc = dvxc + v2rho2
      deallocate (sigma, vsigma, v2rhosigma, v2sigma2)
    case default
      write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
      stop
    end select
    if (xc_sep) then
      select case (xc_family2)
      case (XC_FAMILY_LDA)
        call xc_f03_lda_exc_vxc(xc_func2, np, nrho, e, vrho)
        call xc_f03_lda_fxc(xc_func2, np, nrho, v2rho2)
        exc = exc + e
        vxc = vxc + vrho
        dexc = dexc + irho*(vrho - e)
        dvxc = dvxc + v2rho2
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        allocate (sigma(np), vsigma(np), v2rhosigma(np), v2sigma2(np))
        sigma = 0.0_dp
        call xc_f03_gga_exc_vxc(xc_func2, np, nrho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func2, np, nrho, sigma, v2rho2, v2rhosigma, v2sigma2)
        exc = exc + e
        vxc = vxc + vrho
        dexc = dexc + irho*(vrho - e)
        dvxc = dvxc + v2rho2
        deallocate (sigma, vsigma, v2rhosigma, v2sigma2)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not LDA nor GGA"
        stop
      end select
    end if
    deallocate (nrho, irho, e, vrho, v2rho2)

    exc = exc*hartree
    vxc = vxc*hartree
    dexc = dexc*hartree*abohr3
    dvxc = dvxc*hartree*abohr3
  end subroutine xc_calc_lda

  subroutine xc_calc_gga(rho, drho, ddrho, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma)
    real(dp), intent(in) :: rho(:), drho(:, :), ddrho(:, :, :)
    real(dp), intent(out) :: exc(:), vxc(:), dexcrho(:), dexcsigma(:), dvxcrho(:), dvxcsigma(:)
    integer :: ndim, i, j
    integer(int64) :: np
    real(dp), allocatable :: nrho(:), irho(:), sigma(:), laplacian(:), crossed(:), e(:), &
      &                      vrho(:), vsigma(:), v2rho2(:), v2rhosigma(:), v2sigma2(:), &
      &                      v3rho3(:), v3rho2sigma(:), v3rhosigma2(:), v3sigma3(:)
    np = size(rho, kind=int64)
    ndim = size(drho, 2)
    exc = 0.0_dp
    vxc = 0.0_dp
    dexcrho = 0.0_dp
    dvxcrho = 0.0_dp
    dexcsigma = 0.0_dp
    dvxcsigma = 0.0_dp

    select case (xc_family1)
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      allocate (nrho(np), irho(np), sigma(np), laplacian(np), crossed(np), e(np), &
        &       vrho(np), vsigma(np), v2rho2(np), v2rhosigma(np), v2sigma2(np), &
        &       v3rho3(np), v3rho2sigma(np), v3rhosigma2(np), v3sigma3(np))
      nrho = abohr3*rho
      irho = 1.0_dp/max(tolerance, nrho)
      sigma = 0.0_dp
      laplacian = 0.0_dp
      crossed = 0.0_dp
      do i = 1, ndim
        sigma = sigma + abohr8*drho(:, i)*drho(:, i)
        laplacian = laplacian + abohr5*ddrho(:, i, i)
        do j = 1, ndim
          crossed = crossed + abohr13*drho(:, i)*ddrho(:, j, i)*drho(:, j)
        end do
      end do
      call xc_f03_gga_exc_vxc(xc_func1, np, nrho, sigma, e, vrho, vsigma)
      call xc_f03_gga_fxc(xc_func1, np, nrho, sigma, v2rho2, v2rhosigma, v2sigma2)
      call xc_f03_gga_kxc(xc_func1, np, nrho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      exc = exc + e
      vxc = vxc + vrho - 2.0_dp*(vsigma*laplacian + v2rhosigma*sigma + 2.0_dp*v2sigma2*crossed)
      dexcrho = dexcrho + irho*(vrho - e)
      dexcsigma = dexcsigma + irho*vsigma
      dvxcrho = dvxcrho + v2rho2 - 2.0_dp*(v2rhosigma*laplacian + v3rho2sigma*sigma + 2.0_dp*v3rhosigma2*crossed)
      dvxcsigma = dvxcrho - v2rhosigma - 2.0_dp*(v2sigma2*laplacian + v3rhosigma2*sigma + 2.0_dp*v3sigma3*crossed)
    case default
      write (stderr, "(a)") "[ERROR]: selected functional is not GGA"
      stop
    end select

    if (xc_sep) then
      select case (xc_family2)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_f03_gga_exc_vxc(xc_func2, np, nrho, sigma, e, vrho, vsigma)
        call xc_f03_gga_fxc(xc_func2, np, nrho, sigma, v2rho2, v2rhosigma, v2sigma2)
        call xc_f03_gga_kxc(xc_func2, np, nrho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
        exc = exc + e
        vxc = vxc + vrho - 2.0_dp*(vsigma*laplacian + v2rhosigma*sigma + 2.0_dp*v2sigma2*crossed)
        dexcrho = dexcrho + irho*(vrho - e)
        dexcsigma = dexcsigma + irho*vsigma
        dvxcrho = dvxcrho + v2rho2 - 2.0_dp*(v2rhosigma*laplacian + v3rho2sigma*sigma + 2.0_dp*v3rhosigma2*crossed)
        dvxcsigma = dvxcrho - v2rhosigma - 2.0_dp*(v2sigma2*laplacian + v3rhosigma2*sigma + 2.0_dp*v3sigma3*crossed)
        deallocate (nrho, irho, sigma, laplacian, crossed, e, &
          &         vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, &
          &         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      case default
        write (stderr, "(a)") "[ERROR]: selected functional is not GGA"
        stop
      end select
    end if

    exc = exc*hartree
    vxc = vxc*hartree
    dexcrho = dexcrho*hartree*abohr3
    dexcsigma = dexcsigma*hartree*abohr3
    dvxcrho = dvxcrho*hartree*abohr3
    dvxcsigma = dvxcsigma*hartree*abohr3
  end subroutine xc_calc_gga
end module xc
