subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
  use M_system
  use M_fdata, only: ind2c, numz2c, z2cmax, splineint_2c
  implicit none
  integer, intent(in) :: interaction, isub, in1, in2, non2c, ioption ! Derivative or not
  real*8, intent(in)  :: xin
  real*8, intent(out) :: yout, dfdx
  real*8, parameter :: tol=1.0e-5
  integer :: i, jxx, nnum
  real*8 xmax, xmin, x, h, hi
  real*8 :: a, b, c, d

  xmin = 0.0
  dfdx = 0.0
  jxx = ind2c(interaction,isub)
  nnum = numz2c(jxx,in1,in2)
  xmax = z2cmax(jxx,in1,in2)
  h = (xmax - xmin)/(nnum - 1)
  hi = 1.0/h

  ! Save rounding errors from below
  if(xin .le. xmin) then
    if((xmin - xin)/xmin > tol) then
      write (*,*) 'INTERPOLATE 1D xin, xmin = ', xin, xmin
      stop
    end if
    yout = splineint_2c(1,non2c,1,jxx,in1,in2)
    if(ioption .eq. 1) dfdx = hi*splineint_2c(2,non2c,1,jxx,in1,in2)
    return
  end if

  ! Save rounding errors from above
  if(xin .ge. xmax) then
    if((xin - xmax)/xmax > tol) then
      write (*,*) 'INTERPOLATE 1D xin, xmax = ', xin, xmax
      stop
    end if
    a = splineint_2c(1,non2c,nnum-1,jxx,in1,in2)
    b = splineint_2c(2,non2c,nnum-1,jxx,in1,in2)
    c = splineint_2c(3,non2c,nnum-1,jxx,in1,in2)
    d = splineint_2c(4,non2c,nnum-1,jxx,in1,in2)
    yout = a + b + c + d
    if(ioption .eq. 1) dfdx = hi*(b + 2.0*c + 3.0*d)
    return
  end if

  x = hi*(xin - xmin)
  i = int(x) + 1
  x = x - real(i - 1)

  a = splineint_2c(1,non2c,i,jxx,in1,in2)
  b = splineint_2c(2,non2c,i,jxx,in1,in2)
  c = splineint_2c(3,non2c,i,jxx,in1,in2)
  d = splineint_2c(4,non2c,i,jxx,in1,in2)
  yout = a + x*(b + x*(c + x*d))
  if(ioption .eq. 1) dfdx = hi*(b + x*(2.0*c + x*3.0*d))
end subroutine interpolate_1d

