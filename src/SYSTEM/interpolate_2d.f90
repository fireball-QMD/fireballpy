subroutine interpolate_2d (xin, yin, iauxforce, nx, ny, hx, hy, xintegral, Q_L, dQ_Ldx, dQ_Ldy)
  use iso_c_binding
  use M_fdata, only: numXmax, numYmax
  implicit none
  integer(c_long), intent (in) :: iauxforce, nx, ny
  real(c_double), intent (in) :: xin, yin, hx, hy
  real(c_double), intent (in), dimension (numXmax, numYmax) :: xintegral
  real(c_double), intent (out) :: Q_L, dQ_Ldx, dQ_Ldy
  real(c_double), parameter :: tol = 1.0d-5
  integer(c_long) :: ix, ixp, ixn, ixn2, iy, iyp, iyn, iyn2
  real(c_double) :: xmin, ymin, xmax, ymax, x, y, hxi, hyi, &
    & fm1m1, fm10, fm11, fm12, f0m1, f00, f01, f02, &
    & f1m1, f10, f11, f12, f2m1, f20, f21, f22, &
    & cm10, cm11, cm12, cm13, c00, c01, c02, c03, &
    & c10, c11, c12, c13, c20, c21, c22, c23, &
    & c0, c1, c2, c3, bm1, b0, b1, b2

  ! Initialize
  xmin = 0.0d0
  ymin = 0.0d0
  xmax = xmin + hx*(nx - 1)
  ymax = ymin + hy*(ny - 1)
  dQ_Ldx = 0.0d0
  dQ_Ldy = 0.0d0
  hxi = 1.0d0/hx
  hyi = 1.0d0/hy

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
    x = real(nx - 1, c_double)
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
    y = real(ny - 1, c_double)
    iy = ny - 1
  else
    y = hyi*(yin - ymin)
    iy = int(y) + 1
  end if

  ! Scaled variables
  x = x - real(ix - 1, c_double)
  y = y - real(iy - 1, c_double)

  ! Precompute indices
  ixp = ix - 1
  ixn = ix + 1
  ixn2 = ix + 2
  iyp = iy - 1
  iyn = iy + 1
  iyn2 = iy + 2

  ! Get the surrounding 16 points
  ! If we are on boundaries interpolate the points linearly
  f00 = xintegral(ix,iy)
  f01 = xintegral(ix,iyn)
  f10 = xintegral(ixn,iy)
  f11 = xintegral(ixn,iyn)
  if((ix .eq. 1) .and. (iy .eq. 1)) then
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = 2.0d0*f20 - f21
    f1m1 = 2.0d0*f10 - f11
    f0m1 = 2.0d0*f00 - f01
    fm1m1 = 3.0d0*f00 - f10 - f01
    fm10 = 2.0d0*f00 - f10
    fm11 = 2.0d0*f01 - f11
    fm12 = 2.0d0*f02 - f12
  else if((ix .eq. nx-1) .and. (iy .eq. ny-1)) then
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = 2.0d0*fm11 - fm10
    f02 = 2.0d0*f01 - f00
    f12 = 2.0d0*f11 - f10
    f22 = 3.0d0*f11 - f01 - f10
    f21 = 2.0d0*f11 - f01
    f20 = 2.0d0*f10 - f00
    f2m1 = 2.0d0*f1m1 - f0m1
  else if((ix .eq. 1) .and. (iy .eq. ny-1)) then
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = 2.0d0*f0m1 - f1m1
    fm10 = 2.0d0*f00 - f10
    fm11 = 2.0d0*f01 - f11
    fm12 = 3.0d0*f01 - f00 - f11
    f02 = 2.0d0*f01 - f00
    f12 = 2.0d0*f11 - f10
    f22 = 2.0d0*f21 - f20
  else if((ix .eq. nx-1) .and. (iy .eq. 1)) then
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = 2.0d0*f12 - f02
    f21 = 2.0d0*f11 - f01
    f20 = 2.0d0*f10 - f00
    f2m1 = 3.0d0*f10 - f00 - f11
    f1m1 = 2.0d0*f10 - f11
    f0m1 = 2.0d0*f00 - f01
    fm1m1 = 2.0d0*fm10 - fm11
  else if(ix .eq. 1) then
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = 2.0d0*f0m1 - f1m1
    fm10 = 2.0d0*f00 - f10
    fm11 = 2.0d0*f01 - f11
    fm12 = 2.0d0*f02 - f12
  else if(iy .eq. 1) then
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = 2.0d0*f20 - f21
    f1m1 = 2.0d0*f10 - f11
    f0m1 = 2.0d0*f00 - f01
    fm1m1 = 2.0d0*fm10 - fm11
  else if(ix .eq. nx-1) then
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = 2.0d0*f12 - f02
    f21 = 2.0d0*f11 - f01
    f20 = 2.0d0*f10 - f00
    f2m1 = 2.0d0*f1m1 - f0m1
  else if(iy .eq. ny-1) then
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = 2.0d0*fm11 - fm10
    f02 = 2.0d0*f01 - f00
    f12 = 2.0d0*f11 - f10
    f22 = 2.0d0*f21 - f20
  else
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
  end if

  ! Interpolate x
  cm10 = f0m1
  cm11 = -0.5d0*fm1m1 + 0.5d0*f1m1
  cm12 = fm1m1 - 2.5d0*f0m1 + 2.0d0*f1m1 - 0.5d0*f2m1
  cm13 = -0.5d0*fm1m1 + 1.5d0*f0m1 - 1.5d0*f1m1 + 0.5d0*f2m1
  c00 = f00
  c01 = -0.5d0*fm10 + 0.5d0*f10
  c02 = fm10 - 2.5d0*f00 + 2.0d0*f10 - 0.5d0*f20
  c03 = -0.5d0*fm10 + 1.5d0*f00 - 1.5d0*f10 + 0.5d0*f20
  c10 = f01
  c11 = -0.5d0*fm11 + 0.5d0*f11
  c12 = fm11 - 2.5d0*f01 + 2.0d0*f11 - 0.5d0*f21
  c13 = -0.5d0*fm11 + 1.5d0*f01 - 1.5d0*f11 + 0.5d0*f21
  c20 = f02
  c21 = -0.5d0*fm12 + 0.5d0*f12
  c22 = fm12 - 2.5d0*f02 + 2.0d0*f12 - 0.5d0*f22
  c23 = -0.5d0*fm12 + 1.5d0*f02 - 1.5d0*f12 + 0.5d0*f22
  bm1 = cm10 + x*(cm11 + x*(cm12 + x*cm13))
  b0 = c00 + x*(c01 + x*(c02 + x*c03))
  b1 = c10 + x*(c11 + x*(c12 + x*c13))
  b2 = c20 + x*(c21 + x*(c22 + x*c23))

  ! Interpolate y
  c0 = b0
  c1 = -0.5d0*bm1 + 0.5d0*b1
  c2 = bm1 - 2.5d0*b0 + 2.0d0*b1 - 0.5d0*b2
  c3 = -0.5d0*bm1 + 1.5d0*b0 - 1.5d0*b1 + 0.5d0*b2
  Q_L = c0 + y*(c1 + y*(c2 + y*c3))
  if(iauxforce .eq. 1) then
    dQ_Ldy = hyi*(c1 + y*(2.0d0*c2 + y*3.0d0*c3))
    bm1 = cm11 + x*(2.0d0*cm12 + x*3.0d0*cm13)
    b0 = c01 + x*(2.0d0*c02 + x*3.0d0*c03)
    b1 = c11 + x*(2.0d0*c12 + x*3.0d0*c13)
    b2 = c21 + x*(2.0d0*c22 + x*3.0d0*c23)
    c0 = b0
    c1 = -0.5d0*bm1 + 0.5d0*b1
    c2 = bm1 - 2.5d0*b0 + 2.0d0*b1 - 0.5d0*b2
    c3 = -0.5d0*bm1 + 1.5d0*b0 - 1.5d0*b1 + 0.5d0*b2
    dQ_Ldx = hxi*(c0 + y*(c1 + y*(c2 + y*c3)))
  end if
end subroutine interpolate_2d
