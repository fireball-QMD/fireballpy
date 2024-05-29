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
    & f0m1, f00, f01, f02, f1m1, f10, f11, f12, fx0m1, fx00, fx01, fx02, fx1m1, fx10, fx11, fx12, &
    & pm1, p0, p1, p2, pxm1, px0, px1, px2, pxy0, pxy1, py0, py1, &
    & cxm10, cxm11, cxm12, cxm13, cx00, cx01, cx02, cx03, cx10, cx11, cx12, cx13, cx20, cx21, cx22, cx23, &
    & cy0, cy1, cy2, cy3

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

  ! Scaled variables
  x = x - real(ix - 1)
  y = y - real(iy - 1)

  ! Precompute indices
  ixp = ix - 1
  ixn = ix + 1
  ixn2 = ix + 2
  iyp = iy - 1
  iyn = iy + 1
  iyn2 = iy + 2

  f00 = xintegral(ix,iy)
  f01 = xintegral(ix,iyn)
  f10 = xintegral(ixn,iy)
  f11 = xintegral(ixn,iyn)

  ! Compute x derivatives
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

  ! Interpolate in x
  cx00 = f00
  cx01 = fx00
  cx02 = -3.0*f00 + 3.0*f10 - 2.0*fx00 - fx10
  cx03 = 2.0*f00 - 2.0*f10 + fx00 + fx10
  cx10 = f01
  cx11 = fx01
  cx12 = -3.0*f01 + 3.0*f11 - 2.0*fx01 - fx11
  cx13 = 2.0*f01 - 2.0*f11 + fx01 + fx11
  p0 = cx00 + x*(cx01 + x*(cx02 + x*cx03))
  p1 = cx10 + x*(cx11 + x*(cx12 + x*cx13))
  if(iauxforce .eq. 1) then
    px0 = cx01 + x*(2.0*cx02 + 3.0*x*cx03)
    px1 = cx11 + x*(2.0*cx12 + 3.0*x*cx13)
  end if

  ! Compute y derivatives
  if(iy .eq. 1) then
    f02 = xintegral(ix,iyn2)
    f12 = xintegral(ixn,iyn2)
    if(ix .eq. 1) then
      fx02 = xintegral(ixn,iyn2) - xintegral(ix,iyn2)
      fx12 = 0.5*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2))
    else if(ix .eq. nx-1) then
      fx02 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2))
      fx12 = xintegral(ixn,iyn2) - xintegral(ix,iyn2)
    else
      fx02 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2))
      fx12 = 0.5*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2))
    end if
    cx20 = f02
    cx21 = fx02
    cx22 = -3.0*f02 + 3.0*f12 - 2.0*fx02 - fx12
    cx23 = 2.0*f02 - 2.0*f12 + fx02 + fx12
    p2 = cx20 + x*(cx21 + x*(cx22 + x*cx23))
    py0 = p1 - p0
    py1 = 0.5*(p2 - p0)
    if(iauxforce .eq. 1) then
      px2 = cx21 + x*(2.0*cx22 + x*3.0*cx23)
      pxy0 = px1 - px0
      pxy1 = 0.5*(px2 - px0)
    end if
  else if(iy .eq. ny-1) then
    f0m1 = xintegral(ix,iyp)
    f1m1 = xintegral(ixn,iyp)
    if(ix .eq. 1) then
      fx0m1 = xintegral(ixn,iyp) - xintegral(ix,iyp)
      fx1m1 = 0.5*(xintegral(ixn2,iyp) - xintegral(ix,iyp))
    else if(ix .eq. nx-1) then
      fx0m1 = 0.5*(xintegral(ixn,iyp) - xintegral(ixp,iyp))
      fx1m1 = xintegral(ixn,iyp) - xintegral(ix,iyp)
    else
      fx0m1 = 0.5*(xintegral(ixn,iyp) - xintegral(ixp,iyp))
      fx1m1 = 0.5*(xintegral(ixn2,iyp) - xintegral(ix,iyp))
    end if
    cxm10 = f0m1
    cxm11 = fx0m1
    cxm12 = -3.0*f0m1 + 3.0*f1m1 - 2.0*fx0m1 - fx1m1
    cxm13 = 2.0*f0m1 - 2.0*f1m1 + fx0m1 + fx1m1
    pm1 = cxm10 + x*(cxm11 + x*(cxm12 + x*cxm13))
    py0 = 0.5*(p1 - pm1)
    py1 = p1 - p0
    if(iauxforce .eq. 1) then
      pxm1 = cxm11 + x*(2.0*cxm12 + x*3.0*cxm13)
      pxy0 = 0.5*(px1 - pxm1)
      pxy1 = px1 - px0
    end if
  else
    f0m1 = xintegral(ix,iyp)
    f1m1 = xintegral(ixn,iyp)
    f02 = xintegral(ix,iyn2)
    f12 = xintegral(ixn,iyn2)
    if(ix .eq. 1) then
      fx0m1 = xintegral(ixn,iyp) - xintegral(ix,iyp)
      fx1m1 = 0.5*(xintegral(ixn2,iyp) - xintegral(ix,iyp))
      fx02 = xintegral(ixn,iyn2) - xintegral(ix,iyn2)
      fx12 = 0.5*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2))
    else if(ix .eq. nx-1) then
      fx0m1 = 0.5*(xintegral(ixn,iyp) - xintegral(ixp,iyp))
      fx1m1 = xintegral(ixn,iyp) - xintegral(ix,iyp)
      fx02 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2))
      fx12 = xintegral(ixn,iyn2) - xintegral(ix,iyn2)
    else
      fx0m1 = 0.5*(xintegral(ixn,iyp) - xintegral(ixp,iyp))
      fx1m1 = 0.5*(xintegral(ixn2,iyp) - xintegral(ix,iyp))
      fx02 = 0.5*(xintegral(ixn,iyn2) - xintegral(ixp,iyn2))
      fx12 = 0.5*(xintegral(ixn2,iyn2) - xintegral(ix,iyn2))
    end if
    cxm10 = f0m1
    cxm11 = fx0m1
    cxm12 = -3.0*f0m1 + 3.0*f1m1 - 2.0*fx0m1 - fx1m1
    cxm13 = 2.0*f0m1 - 2.0*f1m1 + fx0m1 + fx1m1
    cx20 = f02
    cx21 = fx02
    cx22 = -3.0*f02 + 3.0*f12 - 2.0*fx02 - fx12
    cx23 = 2.0*f02 - 2.0*f12 + fx02 + fx12
    pm1 = cxm10 + x*(cxm11 + x*(cxm12 + x*cxm13))
    pm1 = cxm10 + x*(cxm11 + x*(cxm12 + x*cxm13))
    p2 = cx20 + x*(cx21 + x*(cx22 + x*cx23))
    py0 = 0.5*(p1 - pm1)
    py1 = 0.5*(p2 - p0)
    if(iauxforce .eq. 1) then
      pxm1 = cxm11 + x*(2.0*cxm12 + x*3.0*cxm13)
      px2 = cx21 + x*(2.0*cx22 + x*3.0*cx23)
      pxy0 = 0.5*(px1 - pxm1)
      pxy1 = 0.5*(px2 - px0)
    end if
  end if

  ! Interpolate in y
  cy0 = p0
  cy1 = py0
  cy2 = -3.0*p0 + 3.0*p1 - 2.0*py0 - py1
  cy3 = 2.0*p0 - 2.0*p1 + py0 + py1
  Q_L = cy0 + y*(cy1 + y*(cy2 + y*cy3))
  if(iauxforce .eq. 1) then
    dQ_Ldy = hyi*(cy1 + y*(2.0*cy2 + y*3.0*cy3))
    cy0 = px0
    cy1 = pxy0
    cy2 = -3.0*px0 + 3.0*px1 - 2.0*pxy0 - pxy1
    cy3 = 2.0*px0 - 2.0*px1 + pxy0 + pxy1
    dQ_Ldx = hxi*(cy0 + y*(cy1 + y*(cy2 + y*cy3)))
  end if
end subroutine interpolate_2d
