subroutine DtrescentrosS (isorp, maxtype, in1, in2, indna, x, y, cost, rhat, sighat, bcnax, f3naXa, f3naXb, f3naXc)
 use M_system
 use M_fdata
 implicit none
 integer, intent (in) :: in1
 integer, intent (in) :: in2
 integer, intent (in) :: indna
 integer, intent (in) :: isorp
 integer, intent (in) :: maxtype
 real*8, intent (in) :: cost
 real*8, intent (in) :: x
 real*8, intent (in) :: y
 real*8, intent (in), dimension (3) :: rhat
 real*8, intent (in), dimension (3) :: sighat
 real*8, intent (out), dimension (nsh_max, nsh_max) :: bcnax
 real*8, intent (out), dimension (3, nsh_max, nsh_max) :: f3naXa
 real*8, intent (out), dimension (3, nsh_max, nsh_max) :: f3naXb
 real*8, intent (out), dimension (3, nsh_max, nsh_max) :: f3naXc
 integer imu
 integer iME
 integer index
 integer inu
 integer ix
 integer kforce
 integer nl
 integer nx
 integer ny
 real*8 amt
 real*8 argument
 real*8 bmt
 real*8 cost2
 real*8 dQ_Ldx
 real*8 dQ_Ldy
 real*8 Q_L
 real*8 sint
 real*8 xxmax
 real*8 yymax
 real*8 hx
 real*8 hy
 real*8 xinv
 real*8, dimension (0:ithetamax - 1, MES_max) :: bcnalist
 real*8, dimension (0:ithetamax - 1) :: dp
 real*8, dimension (MES_max) :: dphlist
 real*8, dimension (0:ithetamax - 1, MES_max) :: dxbcnalist
 real*8, dimension (MES_max) :: dxhlist
 real*8, dimension (0:ithetamax - 1, MES_max) :: dybcnalist
 real*8, dimension (MES_max) :: dyhlist
 real*8, dimension (3, nsh_max, nsh_max) :: f3naMa
 real*8, dimension (3, nsh_max, nsh_max) :: f3naMb
 real*8, dimension (MES_max) :: hlist
 real*8, dimension (0:ithetamax - 1) :: p
 real*8, dimension (nsh_max, nsh_max) :: temp
 do inu = 1, nssh(in2)
   do imu = 1, nssh(in1)
     bcnax(imu,inu) = 0.0d0
     do ix = 1,3
        f3naMa(:,imu,inu) = 0.0d0
        f3naMb(:,imu,inu) = 0.0d0
        f3naXa(:,imu,inu) = 0.0d0
        f3naXb(:,imu,inu) = 0.0d0
     enddo
  end do
 end do
 kforce = 1
 if (ntheta .gt. ithetamax) then
  write (*,*) ' too many thetas, in trescentros.f'
  stop
 end if
 index = icon3c(in1,in2,indna)
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
 do iME = 1, index_maxS(in1,in2)
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(0,iME) = Q_L
    dxbcnalist(0,iME) = dQ_Ldx
    dybcnalist(0,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(1,iME) = Q_L
    dxbcnalist(1,iME) = dQ_Ldx
    dybcnalist(1,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(2,iME) = Q_L
    dxbcnalist(2,iME) = dQ_Ldx
    dybcnalist(2,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(3,iME) = Q_L
    dxbcnalist(3,iME) = dQ_Ldx
    dybcnalist(3,iME) = dQ_Ldy
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(4,iME) = Q_L
    dxbcnalist(4,iME) = dQ_Ldx
    dybcnalist(4,iME) = dQ_Ldy
 end do
 argument = 1.0d0 - cost**2 + 1.0d-05
 if (argument .ge. 0.0d0) then
  sint = sqrt(argument)
 else
  write (*,*) '  Trescentros: BAD SQRT ******* ERROR '
  write (*,*) '  *** STANDARD FIXUP, SET SINT = 0 '
  argument = 0.0d0
  sint = 0.0d0
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
 do iME = 1, index_maxS(in1,in2)
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
  if (mvalueS(iME,in1,in2) .eq. 1) then
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
 call recover_S (in1, in2, hlist, bcnax)
 do ix = 1, 3
  if (x .gt. 1.0d-03) then
    xinv = 1/x
  else
    xinv = x/(0.001**2) 
  end if
  amt = (sighat(ix) - cost*rhat(ix))*xinv
  do iME = 1, index_maxS(in1,in2)
   hlist(iME) = rhat(ix)*dxhlist(iME) + amt*dphlist(iME)
  end do
  call recover_S (in1, in2, hlist, temp)
  do inu = 1, nssh(in2)
   do imu = 1, nssh(in1)
    f3naMa(ix,imu,inu) = temp(imu,inu)
   end do
  end do
  bmt = (cost*sighat(ix) - rhat(ix))/y
  do iME = 1, index_maxS(in1,in2)
    hlist(iME) = - sighat(ix)*dyhlist(iME) + bmt*dphlist(iME) - hlist(iME)/2.0d0
  end do
  call recover_S (in1, in2, hlist, temp)
  do inu = 1, nssh(in2)
   do imu = 1, nssh(in1)
    f3naMb(ix,imu,inu) = temp(imu,inu)
   end do
  end do
 end do
 do inu = 1, nssh(in2)
    do imu = 1, nssh(in1)
      f3naXa(:,imu,inu) = - f3naMa(:,imu,inu)
      f3naXb(:,imu,inu) = - f3naMb(:,imu,inu)
      f3naXc(:,imu,inu) = - f3naXa(:,imu,inu) - f3naXb(:,imu,inu)
    end do
 end do
 return
end subroutine DtrescentrosS

