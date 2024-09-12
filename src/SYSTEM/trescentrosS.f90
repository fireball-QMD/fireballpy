subroutine trescentrosS ( isorp, in1, in2, indna, x, y, cost, bcnax)
  use iso_c_binding
  use M_system, only: ithetamax
  use M_fdata, only: hx_den3, hy_den3, numx3c_den3, numy3c_den3, x3cmax_den3, y3cmax_den3, &
    & den3S_01, den3S_02, den3S_03, den3S_04, den3S_05, &
    & ntheta, MES_max, icon3c, index_maxS, mvalueS, nsh_max, nssh
  implicit none
  integer(c_long), intent (in) :: in1
  integer(c_long), intent (in) :: in2
  integer(c_long), intent (in) :: indna
  integer(c_long), intent (in) :: isorp
  real(c_double), intent (in) :: cost
  real(c_double), intent (in) :: x
  real(c_double), intent (in) :: y
  real(c_double), intent (out), dimension (nsh_max, nsh_max) :: bcnax
  integer(c_long) imu
  integer(c_long) iME
  integer(c_long) index
  integer(c_long) inu
  integer(c_long) kforce
  integer(c_long) nl
  integer(c_long) nx
  integer(c_long) ny
  real(c_double) argument
  real(c_double) cost2
  real(c_double) dQ_Ldx
  real(c_double) dQ_Ldy
  real(c_double) Q_L
  real(c_double) sint
  real(c_double) xxmax
  real(c_double) yymax
  real(c_double) hx
  real(c_double) hy
  real(c_double), dimension (0:ithetamax - 1, MES_max) :: bcnalist
  real(c_double), dimension (MES_max) :: hlist
  real(c_double), dimension (0:ithetamax - 1) :: p
  do inu = 1, nssh(in2)
    do imu = 1, nssh(in1)
      bcnax(imu,inu) = 0.0d0
    end do
  end do
  kforce = 0
  index = icon3c(in1,in2,indna)
    hx = hx_den3(isorp,index)
    hy = hy_den3(isorp,index)
    nx = numx3c_den3(isorp,index)
    ny = numy3c_den3(isorp,index)
    xxmax = x3cmax_den3(isorp,index)
    yymax = y3cmax_den3(isorp,index)
  do iME = 1, index_maxS(in1,in2)
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_01(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(0,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_02(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(1,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_03(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(2,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_04(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(3,iME) = Q_L
    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3S_05(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
    bcnalist(4,iME) = Q_L
  end do
  cost2 = cost**2
  argument = 1.0d0 - cost2
  if (argument .lt. 0.0d0) argument = 1.0d-5
  sint = sqrt(argument)
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
