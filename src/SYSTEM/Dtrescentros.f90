subroutine Dtrescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bcnax, f3naXa, f3naXb, f3naXc)
  use iso_c_binding
  use M_system
  use M_fdata
  implicit none
  integer(c_long), intent (in) :: in1
  integer(c_long), intent (in) :: in2
  integer(c_long), intent (in) :: indna
  integer(c_long), intent (in) :: interaction
  integer(c_long), intent (in) :: isorp
  real(c_double), intent (in) :: cost
  real(c_double), intent (in) :: x
  real(c_double), intent (in) :: y
  real(c_double), intent (in), dimension (3, 3) :: eps
  real(c_double), intent (in), dimension (3, 3, 3) :: depsA
  real(c_double), intent (in), dimension (3, 3, 3) :: depsB
  real(c_double), intent (in), dimension (3) :: rhat
  real(c_double), intent (in), dimension (3) :: sighat
  real(c_double), intent (out), dimension (numorb_max, numorb_max) :: bcnax
  real(c_double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXa
  real(c_double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXb
  real(c_double), intent (out), dimension (3, numorb_max, numorb_max) :: f3naXc
  integer(c_long) imu
  integer(c_long) iME
  integer(c_long) index
  integer(c_long) inu
  integer(c_long) ix
  integer(c_long) kforce
  integer(c_long) nl
  real(c_double) amt
  real(c_double) argument
  real(c_double) bmt
  real(c_double) cost2
  real(c_double) dQ_Ldx
  real(c_double) dQ_Ldy
  real(c_double) Q_L
  real(c_double) sint
  real(c_double) hx, hy
  integer(c_long) nx, ny
  real(c_double) xxmax, yymax, xinv
  real(c_double), dimension (0:ntheta - 1, ME3c_max) :: bcnalist
  real(c_double), dimension (numorb_max, numorb_max) :: bcnam
  real(c_double), dimension (0:ntheta - 1) :: dp
  real(c_double), dimension (ME3c_max) :: dphlist
  real(c_double), dimension (0:ntheta - 1, ME3c_max) :: dxbcnalist
  real(c_double), dimension (ME3c_max) :: dxhlist
  real(c_double), dimension (0:ntheta - 1, ME3c_max) :: dybcnalist
  real(c_double), dimension (ME3c_max) :: dyhlist
  real(c_double), dimension (3, numorb_max, numorb_max) :: f3naMa
  real(c_double), dimension (3, numorb_max, numorb_max) :: f3naMb
  real(c_double), dimension (ME3c_max) :: hlist
  real(c_double), dimension (0:ntheta - 1) :: p
  real(c_double), dimension (numorb_max, numorb_max) :: temp
  kforce = 1
  index = icon3c(in1,in2,indna)
  if (interaction .eq. 1) then
   hx = hx_bcna(isorp,index)
   hy = hy_bcna(isorp,index)
   nx = numx3c_bcna(isorp,index)
   ny = numy3c_bcna(isorp,index)
   xxmax = x3cmax_bcna(isorp,index)
   yymax = y3cmax_bcna(isorp,index)
   if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
     write (*,*) ' What the heck is going on in Dtrescentros!!! error!!! '
     write (*,*) ' x = ', x, ' Max of data = ', xxmax
     write (*,*) ' y = ', y, ' Max of data = ', yymax
     stop
   end if
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
  ! if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
  !   write (*,*) ' What the heck is going on in Dtrescentros!!! error!!! '
  !   write (*,*) ' x = ', x, ' Max of data = ', xxmax
  !   write (*,*) ' y = ', y, ' Max of data = ', yymax
  !   stop
  ! end if
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
   if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
     write (*,*) ' What the heck is going on in trescentros!!! error!!! '
     write (*,*) ' x = ', x, ' Max of data = ', xxmax
     write (*,*) ' y = ', y, ' Max of data = ', yymax
     stop
   end if
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
  if (ntheta .ne. 5) then
    write(*,*) ' ntheta must be 5, but it is ',ntheta
    stop
  end if
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
    if (sint .eq. 0.0d0) then
        write (*,*) ' Dividing by zero (sint = 0) in Dtrescentros.f '
        stop
    end if
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
