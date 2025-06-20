! We assume we do not go outside the range NEVER
subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: ind2c, numz2c, z2cmax, splineint_2c
  implicit none
  integer, intent(in) :: interaction, isub, in1, in2, non2c, ioption ! Derivative or not
  real(double), intent(in)  :: xin
  real(double), intent(out) :: yout, dfdx
  real(double), parameter :: tol=1.0d-5
  integer :: i, jxx, nnum
  real(double) xmax, xmin, x, h, hi
  real(double) :: a, b, c, d

  xmin = 0.0d0
  dfdx = 0.0d0
  jxx = ind2c(interaction,isub)
  nnum = numz2c(jxx,in1,in2)
  xmax = z2cmax(jxx,in1,in2)
  h = (xmax - xmin)/(nnum - 1)
  hi = 1.0d0/h

  ! Save rounding errors from below
  if(xin .le. xmin) then
    yout = splineint_2c(1,non2c,1,jxx,in1,in2)
    if(ioption .eq. 1) dfdx = hi*splineint_2c(2,non2c,1,jxx,in1,in2)
    return
  end if

  ! Save rounding errors from above
  if(xin .ge. xmax) then
    a = splineint_2c(1,non2c,nnum-1,jxx,in1,in2)
    b = splineint_2c(2,non2c,nnum-1,jxx,in1,in2)
    c = splineint_2c(3,non2c,nnum-1,jxx,in1,in2)
    d = splineint_2c(4,non2c,nnum-1,jxx,in1,in2)
    yout = a + b + c + d
    if(ioption .eq. 1) dfdx = hi*(b + 2.0d0*c + 3.0d0*d)
    return
  end if

  x = hi*(xin - xmin)
  i = int(x) + 1
  x = x - real(i - 1, double)

  a = splineint_2c(1,non2c,i,jxx,in1,in2)
  b = splineint_2c(2,non2c,i,jxx,in1,in2)
  c = splineint_2c(3,non2c,i,jxx,in1,in2)
  d = splineint_2c(4,non2c,i,jxx,in1,in2)
  yout = a + x*(b + x*(c + x*d))
  if(ioption .eq. 1) dfdx = hi*(b + x*(2.0d0*c + x*3.0d0*d))
end subroutine interpolate_1d
