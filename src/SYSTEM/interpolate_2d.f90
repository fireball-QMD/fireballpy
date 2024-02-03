subroutine interpolate_2d (xin, yin, iauxforce, nx, ny, hx, hy, xintegral, Q_L, dQ_Ldx, dQ_Ldy)
  use M_system
  use M_constants
  use M_fdata, only: numXmax,numYmax
  implicit none
  integer, intent (in) :: iauxforce
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  real, intent (in) :: xin
  real, intent (in) :: yin
  real, intent (in) :: hx
  real, intent (in) :: hy
  real, intent (in), dimension (numXmax, numYmax) :: xintegral
  real, intent (out) :: Q_L    ! the contibutions for a matrix element
  real, intent (out) :: dQ_Ldx ! d/dx Q_L (Computed only if iauxforce = 1)
  real, intent (out) :: dQ_Ldy ! d/dy Q_L (Computed only if iauxforce = 1)
  real, parameter :: tiny = 1.0d-5
  real, parameter :: small= 1.0d-4
  integer imidx, imidy
  integer k
  real f0p3, f0p6, f1m1, f1m2, f1m3, f1p3, f1p6, f2p1, flm1
  real gradtest
  real gradx, grady
  real px, py
  real, parameter :: xmin = 0
  real, parameter :: ymin = 0
  real bb0,bb1,bb2,bb3
  real, dimension (-1:2,-1:2) :: fun
  real, dimension (-1:2) :: g, gp
  imidx = int((xin - xmin)/hx) + 1
  if (imidx .lt. 2) then
    imidx = 2
  else if (imidx .gt. nx - 2) then
    imidx = nx - 2
  end if
  imidy = int((yin - ymin)/hy) + 1
  if (imidy .lt. 2) then
    imidy = 2
  else if (imidy .gt. ny - 2) then
    imidy = ny - 2
  end if
  px = xin/hx - (imidx - 1)
  py = yin/hy - (imidy - 1)
  fun(-1,-1) = xintegral(imidx - 1,imidy - 1)
  fun(0,-1) = xintegral(imidx,imidy - 1) 
  fun(1,-1) = xintegral(imidx + 1,imidy - 1)
  fun(2,-1) = xintegral(imidx + 2,imidy - 1)
  fun(-1, 0) = xintegral(imidx - 1,imidy)
  fun(0, 0) = xintegral(imidx,imidy) 
  fun(1, 0) = xintegral(imidx + 1,imidy) 
  fun(2, 0) = xintegral(imidx + 2,imidy)
  fun(-1, 1) = xintegral(imidx - 1,imidy + 1)
  fun(0, 1) = xintegral(imidx,imidy + 1)
  fun(1, 1) = xintegral(imidx + 1,imidy + 1)
  fun(2, 1) = xintegral(imidx + 2,imidy + 1)
  fun(-1, 2) = xintegral(imidx - 1,imidy + 2)
  fun(0, 2) = xintegral(imidx,imidy + 2)
  fun(1, 2) = xintegral(imidx + 1,imidy + 2)
  fun(2, 2) = xintegral(imidx + 2,imidy + 2)
  gradx = (fun(1,0) - fun(0,0))/hx
  grady = (fun(0,1) - fun(0,0))/hy
  gradtest = abs(gradx) + abs(grady)
  if (gradtest .lt. tiny) then
   Q_L = (1.0d0 - px - py)*fun(0,0) + px*fun(1,0) + py*fun(0,1)
   if (iauxforce .eq. 1) then
    dQ_Ldx = gradx
    dQ_Ldy = grady
   end if
  else if (gradtest .lt. small) then
   Q_L = py*(py-1)*0.5d0*fun(0,-1) + px*(px-1)*0.5d0*fun(-1,0)   &
     &  + (1 + px*py - px*px - py*py)*fun(0,0)       &
     &  + px*(px - 2*py + 1)*fun(1,0)*0.5d0          &
     &  + py*(py - 2*px + 1)*fun(0,1)*0.5d0          &
     &  + px*py*fun(1,1)
   if (iauxforce .eq. 1) then
    dQ_Ldx = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*py     &
     &      + (fun(-1,0) + fun(1,0) - 2.0d0*fun(0,0))*px         &
     &      - 0.5d0*(fun(-1,0) - fun(1,0)))/hx
    dQ_Ldy = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*px     &
     &      + (fun(0,-1) + fun(0,1) - 2.0d0*fun(0,0))*py         &
     &      - 0.5d0*(fun(0,-1) - fun(0,1)))/hy
   end if
  else if (D2intMeth .eq. 1) then
   do k = -1, 2
    f1m1 = fun(k,-1)
    f1m2 = 2*f1m1
    f1m3 = 3*f1m1
    f0p3 = 3*fun(k,0)
    f0p6 = 2*f0p3
    f1p3 = 3*fun(k,1)
    f1p6 = 2*f1p3
    f2p1 = fun(k,2)
    bb3 = - f1m1 + f0p3 - f1p3 + f2p1
    bb2 = f1m3 - f0p6 + f1p3
    bb1 = - f1m2 - f0p3 + f1p6 - f2p1
    bb0 = f0p6
    g(k) = ((bb3*py + bb2)*py + bb1)*py + bb0
    if (iauxforce .eq. 1) gp(k) = ((3*bb3*py + 2*bb2)*py + bb1)
   end do
   f1m1 = g(-1)
   f1m2 = 2*f1m1
   f1m3 = 3*f1m1
   f0p3 = 3*g(0)
   f0p6 = 2*f0p3
   f1p3 = 3*g(1)
   f1p6 = 2*f1p3
   f2p1 = g(2) 
   bb3 = - f1m1 + f0p3 - f1p3 + f2p1
   bb2 = f1m3 - f0p6 + f1p3
   bb1 = - f1m2 - f0p3 + f1p6 - f2p1
   bb0 = f0p6
   Q_L = (((bb3*px + bb2)*px + bb1)*px + bb0)/36.0d0
   if (iauxforce .eq. 1) then
     dQ_Ldx = ((3*bb3*px + 2*bb2)*px + bb1)/(36.0d0*hx)
     f1m1 = gp(-1)
     f1m2 = 2*f1m1
     f1m3 = 3*f1m1
     f0p3 = 3*gp(0)
     f0p6 = 2*f0p3
     f1p3 = 3*gp(1)
     f1p6 = 2*f1p3
     f2p1 = gp(2) 
     bb3 = - f1m1 + f0p3 - f1p3 + f2p1
     bb2 = f1m3 - f0p6 + f1p3
     bb1 = - f1m2 - f0p3 + f1p6 - f2p1
     bb0 = f0p6
     dQ_Ldy = (((bb3*px + bb2)*px + bb1)*px + bb0)/(36.0d0*hy)
   end if
  else
    write(*,*) ' Invalid interpolation in interpolate_2d '
    write(*,*) ' Some methods that sucked have been removed '
    stop
  end if
  return
  return
end

