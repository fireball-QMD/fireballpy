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
use :: constants, only: abohr3, abohr4, abohr5, abohr6, abohr8, abohr13, hartree, tolerance, err_len
use :: utils, only: utils_safe_inv, utils_log, UTILS_LOG_INFO, UTILS_LOG_ERROR
use :: xc_f03_lib_m, only: xc_func_set_dens_threshold, xc_fxc_f03_func_init_flags, xc_f03_func_get_info, &
    xc_f03_func_info_get_family, xc_f03_func_end, xc_f03_lda_exc_vxc, xc_f03_lda_fxc, xc_f03_lda_kxc, &
    xc_f03_gga_exc_vxc, xc_f03_gga_fxc, xc_f03_gga_kxc, xc_f03_func_t, xc_f03_func_info_t, &
    XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FLAGS_ON_HOST
implicit none
private
public :: xc_isgga, xc_calc, xc_init, xc_end

interface xc_calc
  procedure :: xc_calc_lda
  procedure :: xc_calc_gga
end interface xc_calc

logical :: xc_gga  ! We force either all LDA or all GGA
integer :: nxc
integer, allocatable :: xc_family(:)
real(dp), allocatable :: xc_weights(:)
type(xc_f03_func_t), allocatable :: xc_funcs(:)

contains

pure logical function xc_isgga()
  xc_isgga = xc_gga
end function xc_isgga

subroutine xc_init(names, weights, verbose)
  character(:), allocatable, intent(in) :: names(:)
  real(dp), intent(in) :: weights(:)
  logical, intent(in), optional :: verbose
  integer :: i, iexc, isgga
  type(xc_f03_func_info_t) :: info
  nxc = size(names)
  if (allocated(xc_family)) deallocate(xc_family)
  if (allocated(xc_weights)) deallocate(xc_weights)
  if (allocated(xc_funcs)) deallocate(xc_funcs)
  allocate(xc_family(nxc), xc_weights(nxc), xc_funcs(nxc), xc_info(nxc), isgga(nxc), ws(nxc))

  if (size(weights) /= nxc) then
    call utils_log(UTILS_LOG_ERROR, "XC size(weights) /= size(names)")
    error stop 1
  end if

  isgga = 0
  do i = 1, nxc
    iexc = xc_f03_functional_get_number(trim(names(i)))
    call xc_f03_func_init_flags(xc_funcs(i), iexc, XC_UNPOLARIZED, XC_FLAGS_ON_HOST)
    call xc_func_set_dens_threshold(xc_funcs(i), tolerance)
    xc_weights(i) = weights(i)
    xc_info(i) = xc_f03_func_get_info(xc_funcs(i))
    xc_family(i) = xc_f03_func_info_get_family(xc_info(i))
    select case (xc_family(i))
    case (XC_FAMILY_LDA)
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      isgga = isgga + 1
    case default
      call utils_log(UTILS_LOG_ERROR, "functional "//trim(names(i))//" is not LDA/GGA")
      error stop 1
    end select
  end do
  xc_gga = isgga > 1
end subroutine xc_init

subroutine xc_end()
  integer :: i
  do i = 1, nxc
    call xc_f03_func_end(xc_funcs(i))
  end do
  deallocate(xc_family, xc_weights, xc_funcs)
end subroutine xc_end

subroutine _xc_calc(dens, lapl, grad, hess, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross, &
		   dexcrho2, dvxcrho2)
  real(dp), intent(in) :: dens, lapl
  real(dp), intent(in) :: grad(:), hess(:, :)
  real(dp), intent(out) :: exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross
  real(dp), intent(out), optional :: d2excrho2, d2vxcrho2
  logical :: need2
  integer :: ndim, i, j
  real(dp) :: rho(1), irho(1), sigma(1), laplacian(1), crossed(1), e(1), &
              vrho(1), vsigma(1), v2rho2(1), v2rhosigma(1), v2sigma2(1), &
	      v3rho3(1), v3rho2sigma(1), v3rhosigma2(1), v3sigma3(1)
  need2 = present(dexcrho2) .and. present(dvxcrho2)  ! Will only happen in LDA
  exc = 0.0_dp
  vxc = 0.0_dp
  dexcrho = 0.0_dp
  dvxcrho = 0.0_dp
  dexcsigma = 0.0_dp
  dvxcsigma = 0.0_dp
  dvxclapl = 0.0_dp
  dvxccross = 0.0_dp
  if (need2) then
    d2excrho2 = 0.0_dp
    d2vxcrho2 = 0.0_dp
  end if
  ndim = size(grad)
  rho(1) = abohr3*dens
  irho(1) = utils_safe_inv(rho(1), tolerance)
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
    case (XC_FAMILY_LDA)
      call xc_f03_lda_exc_vxc(xc_funcs(i), 1_int64, rho, e, vrho)
      call xc_f03_lda_fxc(xc_funcs(i), 1_int64, rho, v2rho2)
      if (need2) call xc_f03_lda_kxc(xc_funcs(i), 1_int64, rho, v3rho3)
      exc = exc + xc_weights(i)*e(1)
      vxc = vxc + xc_weights(i)*vrho(1)
      dexcrho = dexcrho + xc_weights(i)*irho(1)*(vrho(1) - e(1))
      dvxcrho = dvxcrho + xc_weights(i)*v2rho2(1)
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call xc_f03_gga_exc_vxc(xc_funcs(i), 1_int64, rho, sigma, e, vrho, vsigma)
      call xc_f03_gga_fxc(xc_funcs(i), 1_int64, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      call xc_f03_gga_kxc(xc_funcs(i), 1_int64, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      exc = exc + xc_weights(i)*e(1)
      vxc = vxc + xc_weights(i)*(vrho(1) - 2.0_dp*(vsigma(1)*laplacian(1) + &
						   v2rhosigma(1)*sigma(1) + &
						   2.0_dp*v2sigma2(1)*crossed(1)))
      dexcrho = dexcrho + xc_weights(i)*irho(1)*(vrho(1) - e(1))
      dexcsigma = dexcsigma + xc_weights(i)*irho(1)*vsigma(1)
      dvxcrho = dvxcrho + xc_weights(i)*(v2rho2(1) - 2.0_dp*(v2rhosigma(1)*laplacian(1) + &
							     v3rho2sigma(1)*sigma(1) + &
							     2.0_dp*v3rhosigma2(1)*crossed(1)))
      dvxcsigma = dvxcsigma - xc_weights(i)*(v2rhosigma(1) - 2.0_dp*(v2sigma2(1)*laplacian(1) + &
								     v3rhosigma2(1)*sigma(1) + &
								     2.0_dp*v3sigma3(1)*crossed(1)))
      dvxclapl = dvxclapl - xc_weights(i)*2.0_dp*vsigma(1)
      dvxccross = dvxccross - xc_weights(i)*4.0_dp*v2sigma2(1)
    case default
      error stop 1
    end select
    if (need2) d2vxcrho2 = d2vxcrho2 + xc_weights(i)*v3rho3(1)
  end do
  exc = exc*hartree
  vxc = vxc*hartree
  dexcrho = dexcrho*hartree*abohr3
  dexcsigma = dexcsigma*hartree*abohr8
  dvxcrho = dvxcrho*hartree*abohr3
  dvxcsigma = dvxcsigma*hartree*abohr8
  dvxclapl = dvxclapl*hartree*abohr5
  dvxccross = dvxccross*hartree*abohr13
  if (need2) then
    d2excrho2 = irho(1)*(dvxcrho - 2.0_dp*dexcrho)*abohr3
    d2vxcrho2 = d2vxcrho2*hartree*abohr6
  end if
end subroutine _xc_calc

subroutine xc_calc_lda(dens, exc, vxc, dexc, dvxc, d2exc, d2vxc)
  real(dp), intent(in) :: dens
  real(dp), intent(out) :: exc, vxc, dexc, dvxc
  real(dp), intent(out), optional :: d2exc, d2vxc
  real(dp) :: lapl, dexcsigma, dvxcsigma, dvxclapl, dvxccross
  real(dp) :: grad(1), hess(1, 1)
  lapl = 0.0_dp
  grad = 0.0_dp
  hess = 0.0_dp
  call _xc_calc(dens, lapl, grad, hess, exc, vxc, dexc, dexcsigma, dvxc, dvxcsigma, dvxclapl, dvxccross, d2exc, d2vxc)
end subroutine xc_calc_lda

subroutine xc_calc_gga(dens, lapl, grad, hess, exc, vxc, dexc, dexcsigma, dvxc, dvxcsigma, dvxclapl, dvxccross)
  real(dp), intent(in) :: dens, lapl
  real(dp), intent(in) :: grad(:), hess(:, :)
  real(dp), intent(out) :: exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross
  call _xc_calc(dens, lapl, grad, hess, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross)
end subroutine xc_calc_gga

end module xc
