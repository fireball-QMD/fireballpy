subroutine smoother (r, rend, xi, stn, dstn)
  use, intrinsic :: iso_fortran_env, only: double => real64
  implicit none
  real(double), intent (in) :: r
  real(double), intent (in) :: rend
  real(double), intent (in) :: xi
  real(double), intent (out) :: stn
  real(double), intent (out) :: dstn
  integer, parameter :: npower = 2
  integer, parameter :: mpower = 2
  integer, parameter :: scaler = 0
  logical, parameter :: old_method = .true.
  real(double) frac
  real(double) rbegin
  real(double) dum
  rbegin = xi*rend
  if (r .gt. rend) then
    stn = 0.0d0
    dstn = 0.0d0
  else if (r .lt. rbegin) then
    stn = 1.0d0
    dstn = 0.0d0
  else
    frac = (r - rbegin)/(rend - rbegin)
    if (old_method) then
      dum = 1.0d0 - frac**npower
      stn = dum**mpower
      dstn = - (mpower*npower)*(dum**(mpower - 1))*(frac**(npower - 1))/(rend - rbegin)
    else  ! new method
      stn = 1 + (scaler-3)*frac**2 + (2-2*scaler)*frac**3 + scaler*frac**4
      dstn=   2*(scaler-3)*frac  + 3*(2-2*scaler)*frac**2 + scaler*frac**3
    end if
  end if
  return
end subroutine smoother
