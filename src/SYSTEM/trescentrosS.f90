subroutine trescentrosS ( isorp, maxtype, in1, in2, indna, x, y, cost, eps, bcnax)
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
  real*8, intent (in), dimension (3, 3) :: eps
  real*8, intent (out), dimension (nsh_max, nsh_max) :: bcnax
  integer imu
  integer iME
  integer index
  integer inu
  integer kforce
  integer nl
  integer nx
  integer ny
  real*8 argument
  real*8 cost2
  real*8 dQ_Ldx
  real*8 dQ_Ldy
  real*8 Q_L
  real*8 sint
  real*8 xxmax
  real*8 yymax
  real*8 hx
  real*8 hy
  real*8, dimension (0:ithetamax - 1, MES_max) :: bcnalist
  real*8, dimension (MES_max) :: hlist
  real*8, dimension (0:ithetamax - 1) :: p
  do inu = 1, nssh(in2)
    do imu = 1, nssh(in1)
      bcnax(imu,inu) = 0.0d0
    end do
  end do
  kforce = 0
  if (ntheta .gt. ithetamax .or. ntheta .gt. 5) then
    write (*,*) ' too many thetas in trescentros.f90 '
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
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(1,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(2,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(3,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(4,iME) = Q_L
  end do
  cost2 = cost**2
  argument = 1.0d0 - cost2 + 1.0d-05
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
  p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
  p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
  p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
  do iME = 1, index_maxS(in1,in2)
    hlist(iME) = 0.0d0
    do nl = 0, ntheta - 1
      hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
    end do
    if (mvalueS(iME,in1,in2) .eq. 1) then
      hlist(iME) = hlist(iME)*sint
    end if
  end do
  call recover_S (in1, in2, hlist, bcnax)
  return
end subroutine trescentrosS

