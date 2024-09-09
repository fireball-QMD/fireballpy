subroutine trescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, bcnax)
  use iso_c_binding
  use M_system, only: numorb_max
  use M_fdata, only: hx_bcna, hy_bcna, hx_den3, hy_den3, numx3c_bcna, numy3c_bcna, &
    & numx3c_den3, numy3c_den3, x3cmax_bcna, y3cmax_bcna, x3cmax_den3, y3cmax_den3, &
    & bcna_01, bcna_02, bcna_03, bcna_04, bcna_05, &
    & den3_01, den3_02, den3_03, den3_04, den3_05, &
    & ntheta, ME3c_max, icon3c, index_max3c, mvalue
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
  real(c_double), intent (out), dimension (numorb_max, numorb_max) :: bcnax
  integer(c_long) iME
  integer(c_long) index
  integer(c_long) kforce
  integer(c_long) nl
  real(c_double) argument
  real(c_double) cost2
  real(c_double) dQ_Ldx
  real(c_double) dQ_Ldy
  real(c_double) Q_L
  real(c_double) sint
  integer(c_long) nx, ny
  real(c_double) xxmax, yymax
  real(c_double) hx, hy
  real(c_double), dimension (0:ntheta - 1, ME3c_max) :: bcnalist
  real(c_double), dimension (numorb_max, numorb_max) :: bcnam
  real(c_double), dimension (ME3c_max) :: hlist
  real(c_double), dimension (0:ntheta - 1) :: p
  kforce = 0
  index = icon3c(in1,in2,indna)
  if (interaction .eq. 1) then
    hx = hx_bcna(isorp,index)
    hy = hy_bcna(isorp,index)
    nx = numx3c_bcna(isorp,index)
    ny = numy3c_bcna(isorp,index)
    xxmax = x3cmax_bcna(isorp,index)
    yymax = y3cmax_bcna(isorp,index)
    if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
       write (*,*) ' What the heck is going on in trescentros!!! error!!! '
       write (*,*) ' x = ', x, ' Max of data = ', xxmax
       write (*,*) ' y = ', y, ' Max of data = ', yymax
       stop
    end if
    do iME = 1, index_max3c(in1,in2)
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_01(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(0,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_02(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(1,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_03(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(2,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_04(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(3,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_05(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(4,iME) = Q_L
    end do
  !else if (interaction .eq. 2) then
  !  hx = hx_xc3c(isorp,index)
  !  hy = hy_xc3c(isorp,index)
  !  nx = numx3c_xc3c(isorp,index)
  !  ny = numy3c_xc3c(isorp,index)
  !  xxmax = x3cmax_xc3c(isorp,index)
  !  yymax = y3cmax_xc3c(isorp,index)
  !  if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
  !     write (*,*) ' What the heck is going on in trescentros!!! error!!! ' 
  !     write (*,*) ' x = ', x, ' Max of data = ', xxmax
  !     write (*,*) ' y = ', y, ' Max of data = ', yymax
  !     stop
  !  end if
  !  do iME = 1, index_max3c(in1,in2)
  !    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !    bcnalist(0,iME) = Q_L
  !    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !    bcnalist(1,iME) = Q_L
  !    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !    bcnalist(2,iME) = Q_L
  !    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !    bcnalist(3,iME) = Q_L
  !    call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
  !    bcnalist(4,iME) = Q_L
  !  end do
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
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_01(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(0,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_02(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(1,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_03(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(2,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_04(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(3,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_05(:,:,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(4,iME) = Q_L
    end do
  end if
  cost2 = cost**2
  argument = 1.0d0 - cost2
  if (argument .lt. 1.0d-5) argument = 1.0d-5
  sint = sqrt(argument) 
  p(0) = 1.0d0
  p(1) = cost
  p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
  p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
  p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
  do iME = 1, index_max3c(in1,in2)
    hlist(iME) = 0.0d0
    do nl = 0, ntheta - 1
      hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
    end do
    if (mvalue(iME,in1,in2) .eq. 1) then
      hlist(iME) = hlist(iME)*sint
    end if
  end do
  call recover_3c (in1, in2, hlist, bcnam)
  call rotate_fb (in1, in2, eps, bcnam, bcnax)
  return
end subroutine trescentros
