subroutine Dtrescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bcnax, f3naXa, f3naXb, f3naXc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max
  use M_fdata, only: hx_bcna, hy_bcna, hx_den3, hy_den3, numx3c_bcna, numy3c_bcna, &
    & numx3c_den3, numy3c_den3, x3cmax_bcna, y3cmax_bcna, x3cmax_den3, y3cmax_den3, &
    & bcna_01, bcna_02, bcna_03, bcna_04, bcna_05, &
    & den3_01, den3_02, den3_03, den3_04, den3_05, &
    & ntheta, ME3c_max, icon3c, index_max3c, mvalue, num_orb
  implicit none
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: indna
  integer, intent (in) :: interaction
  integer, intent (in) :: isorp
  real(double), intent (in) :: cost
  real(double), intent (in) :: x
  real(double), intent (in) :: y
  real(double), intent (in), dimension (3, 3) :: eps
  real(double), intent (in), dimension (3, 3, 3) :: depsA
  real(double), intent (in), dimension (3, 3, 3) :: depsB
  real(double), intent (in), dimension (3) :: rhat
  real(double), intent (in), dimension (3) :: sighat
  real(double), intent (out), dimension (numorb_max, numorb_max) :: bcnax
  real(double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXa
  real(double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXb
  real(double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXc
  integer imu
  integer iME
  integer index
  integer inu
  integer ix
  integer kforce
  integer nl
  real(double) amt
  real(double) argument
  real(double) bmt
  real(double) cost2
  real(double) dQ_Ldx
  real(double) dQ_Ldy
  real(double) Q_L
  real(double) sint
  real(double) hx, hy
  integer nx, ny
  real(double) xxmax, yymax, xinv
  real(double), dimension (0:ntheta - 1, ME3c_max) :: bcnalist
  real(double), dimension (numorb_max, numorb_max) :: bcnam
  real(double), dimension (0:ntheta - 1) :: dp
  real(double), dimension (ME3c_max) :: dphlist
  real(double), dimension (0:ntheta - 1, ME3c_max) :: dxbcnalist
  real(double), dimension (ME3c_max) :: dxhlist
  real(double), dimension (0:ntheta - 1, ME3c_max) :: dybcnalist
  real(double), dimension (ME3c_max) :: dyhlist
  real(double), dimension (3, numorb_max, numorb_max) :: f3naMa
  real(double), dimension (3, numorb_max, numorb_max) :: f3naMb
  real(double), dimension (ME3c_max) :: hlist
  real(double), dimension (0:ntheta - 1) :: p
  real(double), dimension (numorb_max, numorb_max) :: temp
  kforce = 1
  index = icon3c(in1,in2,indna)
  if (interaction .eq. 1) then
   hx = hx_bcna(isorp,index)
   hy = hy_bcna(isorp,index)
   nx = numx3c_bcna(isorp,index)
   ny = numy3c_bcna(isorp,index)
   xxmax = x3cmax_bcna(isorp,index)
   yymax = y3cmax_bcna(isorp,index)
   do iME = 1, index_max3c(in1,in2)
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_01(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(0,iME) = Q_L
    dxbcnalist(0,iME) = dQ_Ldx
    dybcnalist(0,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  bcna_02(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
    bcnalist(1,iME) = Q_L
    dxbcnalist(1,iME) = dQ_Ldx
    dybcnalist(1,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  bcna_03(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
    bcnalist(2,iME) = Q_L
    dxbcnalist(2,iME) = dQ_Ldx
    dybcnalist(2,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  bcna_04(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
    bcnalist(3,iME) = Q_L
    dxbcnalist(3,iME) = dQ_Ldx
    dybcnalist(3,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  bcna_05(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
    bcnalist(4,iME) = Q_L
    dxbcnalist(4,iME) = dQ_Ldx
    dybcnalist(4,iME) = dQ_Ldy
   end do
  !else if (interaction .eq. 2) then
  ! hx = hx_xc3c(isorp,index)
  ! hy = hy_xc3c(isorp,index)
  ! nx = numx3c_xc3c(isorp,index)
  ! ny = numy3c_xc3c(isorp,index)
  ! xxmax = x3cmax_xc3c(isorp,index)
  ! yymax = y3cmax_xc3c(isorp,index)
  ! do iME = 1, index_max3c(in1,in2)
  !  call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  xc3c_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !  bcnalist(0,iME) = Q_L
  !  dxbcnalist(0,iME) = dQ_Ldx
  !  dybcnalist(0,iME) = dQ_Ldy
  !  call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !  bcnalist(1,iME) = Q_L
  !  dxbcnalist(1,iME) = dQ_Ldx
  !  dybcnalist(1,iME) = dQ_Ldy
  !  call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !  bcnalist(2,iME) = Q_L
  !  dxbcnalist(2,iME) = dQ_Ldx
  !  dybcnalist(2,iME) = dQ_Ldy
  !  call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !  bcnalist(3,iME) = Q_L
  !  dxbcnalist(3,iME) = dQ_Ldx
  !  dybcnalist(3,iME) = dQ_Ldy
  !  call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !  bcnalist(4,iME) = Q_L
  !  dxbcnalist(4,iME) = dQ_Ldx
  !  dybcnalist(4,iME) = dQ_Ldy
  ! end do
  else if (interaction .eq. 3) then
   hx = hx_den3(isorp,index)
   hy = hy_den3(isorp,index)
   nx = numx3c_den3(isorp,index)
   ny = numy3c_den3(isorp,index)
   xxmax = x3cmax_den3(isorp,index)
   yymax = y3cmax_den3(isorp,index)
   do iME = 1, index_max3c(in1,in2)
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &       den3_01(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(0,iME) = Q_L
    dxbcnalist(0,iME) = dQ_Ldx
    dybcnalist(0,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_02(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(1,iME) = Q_L
    dxbcnalist(1,iME) = dQ_Ldx
    dybcnalist(1,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  den3_03(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(2,iME) = Q_L
    dxbcnalist(2,iME) = dQ_Ldx
    dybcnalist(2,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_04(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(3,iME) = Q_L
    dxbcnalist(3,iME) = dQ_Ldx
    dybcnalist(3,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_05(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(4,iME) = Q_L
    dxbcnalist(4,iME) = dQ_Ldx
    dybcnalist(4,iME) = dQ_Ldy
   end do
  end if
  argument = 1.0d0 - cost**2
  if (argument .lt. 1.0d-5) argument = 1.0d-5
  sint = sqrt(argument) 
  p(0) = 1.0d0
  p(1) = cost
  cost2 = cost*cost
  p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
  p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
  p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
  dp(0) = 0.0d0
  dp(1) = 1.0d0
  dp(2) = 3.0d0*cost
  dp(3) = (15.0d0*cost2 - 3.0d0)/2.0d0
  dp(4) = (35.0d0*cost*cost2 - 15.0d0*cost)/2.0d0
  do iME = 1, index_max3c(in1,in2)
   hlist(iME) = 0.0d0
   dphlist(iME) = 0.0d0
   dxhlist(iME) = 0.0d0
   dyhlist(iME) = 0.0d0
   do nl = 0, ntheta - 1
    hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
    dphlist(iME) = dphlist(iME) + dp(nl)*bcnalist(nl,iME)
    dxhlist(iME) = dxhlist(iME) + p(nl)*dxbcnalist(nl,iME)
    dyhlist(iME) = dyhlist(iME) + p(nl)*dybcnalist(nl,iME)
   end do
   if (mvalue(iME,in1,in2) .eq. 1) then
    dxhlist(iME) = dxhlist(iME)*sint
    dyhlist(iME) = dyhlist(iME)*sint
    dphlist(iME) = dphlist(iME)*sint - cost*hlist(iME)/sint
    hlist(iME) = hlist(iME)*sint
   end if
  end do
  call recover_3c (in1, in2, hlist, bcnam)
  call rotate_fb (in1, in2,  eps, bcnam, bcnax)
  do ix = 1, 3
   if (x .gt. 1.0d-03) then
     xinv = 1/x
   else
     xinv = x/(0.001**2) 
   end if
   amt = (sighat(ix) - cost*rhat(ix))*xinv
   do iME = 1, index_max3c(in1,in2)
    hlist(iME) = rhat(ix)*dxhlist(iME) + amt*dphlist(iME)
   end do
   call recover_3c (in1, in2, hlist, temp)
   do inu = 1, num_orb(in2)
    do imu = 1, num_orb(in1)
     f3naMa(ix,imu,inu) = temp(imu,inu)
    end do
   end do
   bmt = (cost*sighat(ix) - rhat(ix))/y
   do iME = 1, index_max3c(in1,in2)
     hlist(iME) = - sighat(ix)*dyhlist(iME) + bmt*dphlist(iME) - hlist(iME)/2.0d0
   end do
   call recover_3c (in1, in2, hlist, temp)
   do inu = 1, num_orb(in2)
    do imu = 1, num_orb(in1)
     f3naMb(ix,imu,inu) = temp(imu,inu)
    end do
   end do
  end do
  call rotated (in1, in2, eps, depsA, bcnam, f3naMa, f3naXa)
  call rotated (in1, in2, eps, depsB, bcnam, f3naMb, f3naXb)
  do inu = 1, num_orb(in2)
    do imu = 1, num_orb(in1)
     f3naXa(:,imu,inu) = - f3naXa(:,imu,inu)
     f3naXb(:,imu,inu) = - f3naXb(:,imu,inu)
     f3naXc(:,imu,inu) = - f3naXa(:,imu,inu) - f3naXb(:,imu,inu)
    end do
  end do
  return
end subroutine Dtrescentros
