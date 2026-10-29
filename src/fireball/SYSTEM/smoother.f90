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
  integer, parameter :: scaler = 0
  logical, parameter :: old_method = .true.
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
    if (old_method) then
      dum = 1.0_dp - frac**npower
      stn = dum**mpower
      dstn = - (mpower*npower)*(dum**(mpower - 1))*(frac**(npower - 1))/(rend - rbegin)
    else  ! new method
      stn = 1 + (scaler-3)*frac**2 + (2-2*scaler)*frac**3 + scaler*frac**4
      dstn=   2*(scaler-3)*frac  + 3*(2-2*scaler)*frac**2 + scaler*frac**3
    end if
  end if
  return
end subroutine smoother
