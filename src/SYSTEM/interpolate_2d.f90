subroutine interpolate_2d (xin, yin, iauxforce, nx, ny, hx, hy, xintegral, Q_L, dQ_Ldx, dQ_Ldy)
  use M_fdata, only: numXmax, numYmax
  implicit none
  integer, intent (in) :: iauxforce, nx, ny
  real*8, intent (in) :: xin, yin, hx, hy
  real*8, intent (in), dimension (numXmax, numYmax) :: xintegral
  real*8, intent (out) :: Q_L, dQ_Ldx, dQ_Ldy
  real*8, parameter :: tol = 1.0e-5
  integer :: ix, ixp, ixn, ixn2, iy, iyp, iyn, iyn2
  real*8 :: xmin, ymin, xmax, ymax, x, y, hxi, hyi, &
    & f00, f01, f10, f11, fx00, fx01, fx10, fx11, &
    & fy00, fy01, fy10, fy11, fxy00, fxy01, fxy10, fxy11
  real*8, dimension(4,4) :: coefs

  ! Initialize
  xmin = 0.0
  ymin = 0.0
  xmax = xmin + hx*(nx - 1)
  ymax = ymin + hy*(ny - 1)
  dQ_Ldx = 0.0
  dQ_Ldy = 0.0
  hxi = 1.0/hx
  hyi = 1.0/hy

  ! Save rounding errors for the x
  if(xin .le. xmin) then
    if((xmin - xin)/xmin > tol) then
      write (*,*) 'INTERPOLATE 2D xin, xmin = ', xin, xmin
      stop
    end if
    x = 0.0
    ix = 1
  else if(xin .ge. xmax) then
    if((xin - xmax)/xmax > tol) then
      write (*,*) 'INTERPOLATE 2D xin, xmax = ', xin, xmax
      stop
    end if
    x = real(nx - 1)
    ix = nx - 1
  else
    x = hxi*(xin - xmin)
    ix = int(x) + 1
  end if

  ! Save rounding errors for the y
  if(yin .le. ymin) then
    if((ymin - yin)/ymin > tol) then
      write (*,*) 'INTERPOLATE 2D yin, ymin = ', yin, ymin
      stop
    end if
    y = 0.0
    iy = 1
  else if(yin .ge. ymax) then
    if((yin - ymax)/ymax > tol) then
      write (*,*) 'INTERPOLATE 2D yin, ymax = ', yin, ymax
      stop
    end if
    y = real(ny - 1)
    iy = ny - 1
  else
    y = hyi*(yin - ymin)
    iy = int(y) + 1
  end if

  ! Precompute indices
  ixp = ix - 1
  ixn = ix + 1
  ixn2 = ix + 2
  iyp = iy - 1
  iyn = iy + 1
  iyn2 = iy + 2

  ! Get necessary coefficients with care of the boundaries
  ! f(x,y)
  f00 = xintegral(ix,iy)
  f01 = xintegral(ix,iyn)
  f10 = xintegral(ixn,iy)
  f11 = xintegral(ixn,iyn)
  ! df(x,y)/dx
  if(ix .eq. 1) then
    fx00 = xintegral(ixn,iy) - xintegral(ix,iy)
    fx01 = xintegral(ixn,iyn) - xintegral(ix,iyn)
    fx10 = 0.5*(xintegral(ixn2,iy) - xintegral(ix,iy))
    fx11 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn))
  else if(ix .eq. nx-1) then
    fx00 = 0.5*(xintegral(ixn,iy) - xintegral(ixp,iy))
    fx01 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn))
    fx10 = xintegral(ixn,iy) - xintegral(ix,iy)
    fx11 = xintegral(ixn,iyn) - xintegral(ix,iyn)
  else
    fx00 = 0.5*(xintegral(ixn,iy) - xintegral(ixp,iy))
    fx01 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn))
    fx10 = 0.5*(xintegral(ixn2,iy) - xintegral(ix,iy))
    fx11 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn))
  end if
  ! df(x,y)/dy
  if(iy .eq. 1) then
    fy00 = xintegral(ix,iyn) - xintegral(ix,iy)
    fy01 = 0.5*(xintegral(ix,iyn2) - xintegral(ix,iy))
    fy10 = xintegral(ixn,iyn) - xintegral(ixn,iy)
    fy11 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixn,iy))
  else if(iy .eq. ny-1) then
    fy00 = 0.5*(xintegral(ix,iyn) - xintegral(ix,iyp))
    fy01 = xintegral(ix,iyn) - xintegral(ix,iy)
    fy10 = 0.5*(xintegral(ixn,iyn) - xintegral(ixn,iyp))
    fy11 = xintegral(ixn,iyn) - xintegral(ixn,iy)
  else
    fy00 = 0.5*(xintegral(ix,iyn) - xintegral(ix,iyp))
    fy01 = 0.5*(xintegral(ix,iyn2) - xintegral(ix,iy))
    fy10 = 0.5*(xintegral(ixn,iyn) - xintegral(ixn,iyp))
    fy11 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixn,iy))
  end if
  ! d^2f(x,y)/dxdy
  if((ix .eq. 1) .and. (iy .eq. 1)) then
    fxy00 = xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iy) + xintegral(ix,iy)
    fxy01 = 0.5*(xintegral(ixn,iyn2) - xintegral(ix,iyn2) - xintegral(ixn,iy) + xintegral(ix,iy))
    fxy10 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iy) + xintegral(ix,iy))
    fxy11 = 0.25*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2) - xintegral(ixn2,iy) + xintegral(ix,iy))
  else if((ix .eq. nx-1) .and. (iy .eq. ny-1)) then
    fxy00 = 0.25*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iyp) + xintegral(ixp,iyp))
    fxy01 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iy) + xintegral(ixp,iy))
    fxy10 = 0.5*(xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iyp) + xintegral(ix,iyp))
    fxy11 = xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iy) + xintegral(ix,iy)
  else if((ix .eq. 1) .and. (iy .eq. ny-1)) then
    fxy00 = 0.5*(xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iyp) + xintegral(ix,iyp))
    fxy01 = xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iy) + xintegral(ix,iy)
    fxy10 = 0.25*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iyp) + xintegral(ix,iyp))
    fxy11 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iy) + xintegral(ix,iy))
  else if((ix .eq. nx-1) .and. (iy .eq. 1)) then
    fxy00 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iy) + xintegral(ixp,iy))
    fxy01 = 0.25*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2) - xintegral(ixn,iy) + xintegral(ixp,iy))
    fxy10 = xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iy) + xintegral(ix,iy)
    fxy11 = 0.5*(xintegral(ixn,iyn2) - xintegral(ix,iyn2) - xintegral(ixn,iy) + xintegral(ix,iy))
  else if(ix .eq. 1) then
    fxy00 = 0.5*(xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iyp) + xintegral(ix,iyp))
    fxy01 = 0.5*(xintegral(ixn,iyn2) - xintegral(ix,iyn2) - xintegral(ixn,iy)  + xintegral(ix,iy))
    fxy10 = 0.25*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iyp) + xintegral(ix,iyp))
    fxy11 = 0.25*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2) - xintegral(ixn2,iy) + xintegral(ix,iy))
  else if(iy .eq. 1) then
    fxy00 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iy) + xintegral(ixp,iy))
    fxy01 = 0.25*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2) - xintegral(ixn,iy)  + xintegral(ixp,iy))
    fxy10 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iy) + xintegral(ix,iy))
    fxy11 = 0.25*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2) - xintegral(ixn2,iy) + xintegral(ix,iy))
  else if(ix .eq. nx-1) then
    fxy00 = 0.25*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iyp) + xintegral(ixp,iyp))
    fxy01 = 0.25*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2) - xintegral(ixn,iy)  + xintegral(ixp,iy))
    fxy10 = 0.5*(xintegral(ixn,iyn) - xintegral(ix,iyn) - xintegral(ixn,iyp) + xintegral(ix,iyp))
    fxy11 = 0.5*(xintegral(ixn,iyn2) - xintegral(ix,iyn2) - xintegral(ixn,iy) + xintegral(ix,iy))
  else if(iy .eq. ny-1) then
    fxy00 = 0.25*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iyp) + xintegral(ixp,iyp))
    fxy01 = 0.5*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iy)  + xintegral(ixp,iy))
    fxy10 = 0.25*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iyp) + xintegral(ix,iyp))
    fxy11 = 0.5*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iy) + xintegral(ix,iy))
  else
    fxy00 = 0.25*(xintegral(ixn,iyn) - xintegral(ixp,iyn) - xintegral(ixn,iyp) + xintegral(ixp,iyp))
    fxy01 = 0.25*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2) - xintegral(ixn,iy)  + xintegral(ixp,iy))
    fxy10 = 0.25*(xintegral(ixn2,iyn) - xintegral(ix,iyn) - xintegral(ixn2,iyp) + xintegral(ix,iyp))
    fxy11 = 0.25*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2) - xintegral(ixn2,iy) + xintegral(ix,iy))
  end if

  ! Coefficient matrix
  coefs(1,1) = f00
  coefs(1,2) = fy00
  coefs(1,3) = -3.0*f00 + 3.0*f01 - 2.0*fy00 - fy01
  coefs(1,4) = 2.0*f00 - 2.0*f01 + fy00 + fy01
  coefs(2,1) = fx00
  coefs(2,2) = fxy00
  coefs(2,3) = -3.0*fx00 + 3.0*fx01 - 2.0*fxy00 - fxy01
  coefs(2,4) = 2.0*fx00 - 2.0*fx01 + fxy00 + fxy01
  coefs(3,1) = -3.0*f00 + 3.0*f10 - 2.0*fx00 - fx10
  coefs(3,2) = -3.0*fy00 + 3.0*fy10 - 2.0*fxy00 - fxy10
  coefs(3,3) = 9.0*f00 - 9.0*f01 + 6.0*fy00 + 3.0*fy01 - 9.0*f10 + 9.0*f11 - 6.0*fy10 - 3.0*fy11 + 6.0*fx00 &
    & - 6.0*fx01 + 4.0*fxy00 + 2.0*fx01 + 3.0*fx10 - 3.0*fx11 + 2.0*fxy10 + fxy11
  coefs(3,4) = -6.0*f00 + 6.0*f01 - 3.0*fy00 - 3.0*fy01 + 6.0*f10 - 6.0*f11 + 3.0*fy10 + 3.0*fy11 - 4.0*fx00 &
    & + 4.0*fx01 - 2.0*fxy00 - 2.0*fx01 - 2.0*fx10 + 2.0*fx11 - fxy10 - fxy11
  coefs(4,1) = 2.0*f00 - 2.0*f10 + fx00 + fx10
  coefs(4,2) = 2.0*fy00 - 2.0*fy10 + fxy00 + fxy10
  coefs(4,3) = -6.0*f00 + 6.0*f01 - 4.0*fy00 - 2.0*fy01 + 6.0*f10 - 6.0*f11 + 4.0*fy10 + 2.0*fy11 - 3.0*fx00 &
    & + 3.0*fx01 - 2.0*fxy00 - fx01 - 3.0*fx10 + 3.0*fx11 - 2.0*fxy10 - fxy11
  coefs(4,4) = 4.0*f00 - 4.0*f01 + 2.0*fy00 + 2.0*fy01 - 4.0*f10 + 4.0*f11 - 2.0*fy10 - 2.0*fy11 + 2.0*fx00 &
    & - 2.0*fx01 + fxy00 + fx01 + 2.0*fx10 - 2.0*fx11 + fxy10 + fxy11

  ! Scaled variables
  x = x - real(ix - 1)
  y = y - real(iy - 1)

  ! Compute polynomial
  Q_L =  coefs(1,1) + y*(coefs(1,2) + y*(coefs(1,3) + y*coefs(1,4))) + &
    & x*(coefs(2,1) + y*(coefs(2,2) + y*(coefs(2,3) + y*coefs(2,4))) + &
    & x*(coefs(3,1) + y*(coefs(3,2) + y*(coefs(3,3) + y*coefs(3,4))) + &
    & x*(coefs(4,1) + y*(coefs(4,2) + y*(coefs(4,3) + y*coefs(4,4))))))

  ! Compute derivatives
  if(iauxforce .eq. 1) then
    dQ_Ldx = hxi*(coefs(2,1) + y*(coefs(2,2) + y*(coefs(2,3) + y*coefs(2,4))) + &
      &    2.0*x*(coefs(3,1) + y*(coefs(3,2) + y*(coefs(3,3) + y*coefs(3,4))) + &
      &    3.0*x*(coefs(4,1) + y*(coefs(4,2) + y*(coefs(4,3) + y*coefs(4,4))))))
    dQ_Ldy = hyi*(coefs(1,2) + 2.0*y*(coefs(1,3) + 3.0*y*coefs(1,4)) + &
      &        x*(coefs(2,2) + 2.0*y*(coefs(2,3) + 3.0*y*coefs(2,4)) + &
      &        x*(coefs(3,2) + 2.0*y*(coefs(3,3) + 3.0*y*coefs(3,4)) + &
      &        x*(coefs(4,2) + 2.0*y*(coefs(4,3) + 3.0*y*coefs(4,4))))))
  end if
end subroutine interpolate_2d
