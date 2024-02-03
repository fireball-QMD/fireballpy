subroutine smoother (r, rend, xi, stn, dstn)
  implicit none
  real, intent (in) :: r
  real, intent (in) :: rend
  real, intent (in) :: xi
  real, intent (out) :: stn
  real, intent (out) :: dstn
  integer, parameter :: npower = 2
  integer, parameter :: mpower = 2
  integer, parameter :: scaler = 0
  logical, parameter :: old_method = .true.
  real frac
  real rbegin
  real dum
  rbegin = xi*rend
  if (r .lt. 0.0d0) then
    write (*,*) ' r < 0 in smoother *** error! '
    stop
  else if (r .gt. rend) then
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
end

