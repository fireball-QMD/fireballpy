! We assume we do not go outside the range NEVER
subroutine interpolate_2d (xin, yin, iauxforce, nx, ny, hx, hy, xintegral, Q_L, dQ_Ldx, dQ_Ldy)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_fdata, only: numXmax, numYmax
  implicit none
  integer, intent (in) :: iauxforce, nx, ny
  real(dp), intent (in) :: xin, yin, hx, hy
  real(dp), intent (in), dimension (numXmax, numYmax) :: xintegral
  real(dp), intent (out) :: Q_L, dQ_Ldx, dQ_Ldy
  real(dp), parameter :: tol = 1.0e-5_dp
  integer :: ix, ixp, ixn, ixn2, iy, iyp, iyn, iyn2
  real(dp) :: xmin, ymin, xmax, ymax, x, y, hxi, hyi, &
    & fm1m1, fm10, fm11, fm12, f0m1, f00, f01, f02, &
    & f1m1, f10, f11, f12, f2m1, f20, f21, f22, &
    & cm10, cm11, cm12, cm13, c00, c01, c02, c03, &
    & c10, c11, c12, c13, c20, c21, c22, c23, &
    & c0, c1, c2, c3, bm1, b0, b1, b2

  ! Initialize
  xmin = 0.0_dp
  ymin = 0.0_dp
  xmax = xmin + hx*(nx - 1)
  ymax = ymin + hy*(ny - 1)
  dQ_Ldx = 0.0_dp
  dQ_Ldy = 0.0_dp
  hxi = 1.0_dp/hx
  hyi = 1.0_dp/hy

  ! Save rounding errors for the x
  if(xin .le. xmin) then
    x = 0.0
    ix = 1
  else if(xin .ge. xmax) then
    x = real(nx - 1, dp)
    ix = nx - 1
  else
    x = hxi*(xin - xmin)
    ix = int(x) + 1
  end if

  ! Save rounding errors for the y
  if(yin .le. ymin) then
    y = 0.0
    iy = 1
  else if(yin .ge. ymax) then
    y = real(ny - 1, dp)
    iy = ny - 1
  else
    y = hyi*(yin - ymin)
    iy = int(y) + 1
  end if

  ! Scaled variables
  x = x - real(ix - 1, dp)
  y = y - real(iy - 1, dp)

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
    f2m1 = 2.0_dp*f20 - f21
    f1m1 = 2.0_dp*f10 - f11
    f0m1 = 2.0_dp*f00 - f01
    fm1m1 = 3.0_dp*f00 - f10 - f01
    fm10 = 2.0_dp*f00 - f10
    fm11 = 2.0_dp*f01 - f11
    fm12 = 2.0_dp*f02 - f12
  else if((ix .eq. nx-1) .and. (iy .eq. ny-1)) then
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = 2.0_dp*fm11 - fm10
    f02 = 2.0_dp*f01 - f00
    f12 = 2.0_dp*f11 - f10
    f22 = 3.0_dp*f11 - f01 - f10
    f21 = 2.0_dp*f11 - f01
    f20 = 2.0_dp*f10 - f00
    f2m1 = 2.0_dp*f1m1 - f0m1
  else if((ix .eq. 1) .and. (iy .eq. ny-1)) then
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = 2.0_dp*f0m1 - f1m1
    fm10 = 2.0_dp*f00 - f10
    fm11 = 2.0_dp*f01 - f11
    fm12 = 3.0_dp*f01 - f00 - f11
    f02 = 2.0_dp*f01 - f00
    f12 = 2.0_dp*f11 - f10
    f22 = 2.0_dp*f21 - f20
  else if((ix .eq. nx-1) .and. (iy .eq. 1)) then
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = 2.0_dp*f12 - f02
    f21 = 2.0_dp*f11 - f01
    f20 = 2.0_dp*f10 - f00
    f2m1 = 3.0_dp*f10 - f00 - f11
    f1m1 = 2.0_dp*f10 - f11
    f0m1 = 2.0_dp*f00 - f01
    fm1m1 = 2.0_dp*fm10 - fm11
  else if(ix .eq. 1) then
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = 2.0_dp*f0m1 - f1m1
    fm10 = 2.0_dp*f00 - f10
    fm11 = 2.0_dp*f01 - f11
    fm12 = 2.0_dp*f02 - f12
  else if(iy .eq. 1) then
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = xintegral(ixn2, iyn2)
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = 2.0_dp*f20 - f21
    f1m1 = 2.0_dp*f10 - f11
    f0m1 = 2.0_dp*f00 - f01
    fm1m1 = 2.0_dp*fm10 - fm11
  else if(ix .eq. nx-1) then
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = xintegral(ixp, iyn2)
    f02 = xintegral(ix, iyn2)
    f12 = xintegral(ixn, iyn2)
    f22 = 2.0_dp*f12 - f02
    f21 = 2.0_dp*f11 - f01
    f20 = 2.0_dp*f10 - f00
    f2m1 = 2.0_dp*f1m1 - f0m1
  else if(iy .eq. ny-1) then
    f21 = xintegral(ixn2, iyn)
    f20 = xintegral(ixn2, iy)
    f2m1 = xintegral(ixn2, iyp)
    f1m1 = xintegral(ixn, iyp)
    f0m1 = xintegral(ix, iyp)
    fm1m1 = xintegral(ixp, iyp)
    fm10 = xintegral(ixp, iy)
    fm11 = xintegral(ixp, iyn)
    fm12 = 2.0_dp*fm11 - fm10
    f02 = 2.0_dp*f01 - f00
    f12 = 2.0_dp*f11 - f10
    f22 = 2.0_dp*f21 - f20
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
  cm11 = -0.5_dp*fm1m1 + 0.5_dp*f1m1
  cm12 = fm1m1 - 2.5_dp*f0m1 + 2.0_dp*f1m1 - 0.5_dp*f2m1
  cm13 = -0.5_dp*fm1m1 + 1.5_dp*f0m1 - 1.5_dp*f1m1 + 0.5_dp*f2m1
  c00 = f00
  c01 = -0.5_dp*fm10 + 0.5_dp*f10
  c02 = fm10 - 2.5_dp*f00 + 2.0_dp*f10 - 0.5_dp*f20
  c03 = -0.5_dp*fm10 + 1.5_dp*f00 - 1.5_dp*f10 + 0.5_dp*f20
  c10 = f01
  c11 = -0.5_dp*fm11 + 0.5_dp*f11
  c12 = fm11 - 2.5_dp*f01 + 2.0_dp*f11 - 0.5_dp*f21
  c13 = -0.5_dp*fm11 + 1.5_dp*f01 - 1.5_dp*f11 + 0.5_dp*f21
  c20 = f02
  c21 = -0.5_dp*fm12 + 0.5_dp*f12
  c22 = fm12 - 2.5_dp*f02 + 2.0_dp*f12 - 0.5_dp*f22
  c23 = -0.5_dp*fm12 + 1.5_dp*f02 - 1.5_dp*f12 + 0.5_dp*f22
  bm1 = cm10 + x*(cm11 + x*(cm12 + x*cm13))
  b0 = c00 + x*(c01 + x*(c02 + x*c03))
  b1 = c10 + x*(c11 + x*(c12 + x*c13))
  b2 = c20 + x*(c21 + x*(c22 + x*c23))

  ! Interpolate y
  c0 = b0
  c1 = -0.5_dp*bm1 + 0.5_dp*b1
  c2 = bm1 - 2.5_dp*b0 + 2.0_dp*b1 - 0.5_dp*b2
  c3 = -0.5_dp*bm1 + 1.5_dp*b0 - 1.5_dp*b1 + 0.5_dp*b2
  Q_L = c0 + y*(c1 + y*(c2 + y*c3))
  if(iauxforce .eq. 1) then
    dQ_Ldy = hyi*(c1 + y*(2.0_dp*c2 + y*3.0_dp*c3))
    bm1 = cm11 + x*(2.0_dp*cm12 + x*3.0_dp*cm13)
    b0 = c01 + x*(2.0_dp*c02 + x*3.0_dp*c03)
    b1 = c11 + x*(2.0_dp*c12 + x*3.0_dp*c13)
    b2 = c21 + x*(2.0_dp*c22 + x*3.0_dp*c23)
    c0 = b0
    c1 = -0.5_dp*bm1 + 0.5_dp*b1
    c2 = bm1 - 2.5_dp*b0 + 2.0_dp*b1 - 0.5_dp*b2
    c3 = -0.5_dp*bm1 + 1.5_dp*b0 - 1.5_dp*b1 + 0.5_dp*b2
    dQ_Ldx = hxi*(c0 + y*(c1 + y*(c2 + y*c3)))
  end if
end subroutine interpolate_2d
