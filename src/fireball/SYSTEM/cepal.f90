! On-site OLSXC/SNXC correction: the LDA of Ceperley-Alder as parameterised by
! Perdew-Zunger, Phys. Rev. B23, 5048 (1981).
!
! The functional is no longer evaluated here: it is delegated to libXC through
! src/common/xc.f90, using XC_LDA_X (1) + XC_LDA_C_PZ (9), which is the same
! functional. The signature is kept unchanged so that none of the call sites in
! the OLSXC assemblers need to be touched.
subroutine cepal (rh, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use :: xc, only: xc_calc, xc_init
  implicit none
  real(double), intent (in) :: rh
  real(double), intent (out) :: exc, muxc, dexc, d2exc, dmuxc, d2muxc
  real(double), parameter :: delta_rh = 1.0d-6
  logical, save :: iniciado = .false.
  real(double) :: rhx

  ! Lazily initialised; the functional here is always the same, so there is no
  ! teardown. Safe because none of the callers run inside a parallel region.
  if (.not. iniciado) then
    call xc_init(1, 9)
    iniciado = .true.
  end if

  ! Same regularisation as the original implementation.
  rhx = sqrt(rh*rh + delta_rh)
  call xc_calc(rhx, exc, muxc, dexc, dmuxc, d2exc, d2muxc)
end subroutine cepal
