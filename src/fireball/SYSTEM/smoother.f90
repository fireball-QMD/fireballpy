subroutine smoother (r, rend, xi, stn, dstn)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent (in) :: r
  real(dp), intent (in) :: rend
  real(dp), intent (in) :: xi
  real(dp), intent (out) :: stn
  real(dp), intent (out) :: dstn
  integer, parameter :: npower = 2
  integer, parameter :: mpower = 2
  real(dp) frac
  real(dp) rbegin
  real(dp) dum
  rbegin = xi*rend
  if (r .gt. rend) then
    stn = 0.0_dp
    dstn = 0.0_dp
  else if (r .lt. rbegin) then
    stn = 1.0_dp
    dstn = 0.0_dp
  else
    frac = (r - rbegin)/(rend - rbegin)
    dum = 1.0_dp - frac**npower
    stn = dum**mpower
    dstn = - (mpower*npower)*(dum**(mpower - 1))*(frac**(npower - 1))/(rend - rbegin)
  end if
  return
end subroutine smoother
