subroutine trescentros (interaction, isorp, maxtype, in1, in2, indna, x, y, cost, eps, bcnax)
  use M_system
  use M_fdata
  implicit none
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: indna
  integer, intent (in) :: interaction
  integer, intent (in) :: isorp
  integer, intent (in) :: maxtype
  real(8), intent (in) :: cost
  real(8), intent (in) :: x
  real(8), intent (in) :: y
  real(8), intent (in), dimension (3, 3) :: eps
  real(8), intent (out), dimension (numorb_max, numorb_max) :: bcnax
  integer imu
  integer iME
  integer index
  integer inu
  integer kforce
  integer nl
  real(8) argument
  real(8) cost2
  real(8) dQ_Ldx
  real(8) dQ_Ldy
  real(8) Q_L
  real(8) sint
  integer nx, ny
  real(8) xxmax, yymax
  real(8) hx, hy
  real(8), dimension (0:ntheta - 1, ME3c_max) :: bcnalist
  real(8), dimension (numorb_max, numorb_max) :: bcnam
  real(8), dimension (ME3c_max) :: hlist
  real(8), dimension (0:ntheta - 1) :: p
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
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(0,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(1,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(2,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(3,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, bcna_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(4,iME) = Q_L
    end do
  else if (interaction .eq. 2) then
    hx = hx_xc3c(isorp,index)
    hy = hy_xc3c(isorp,index)
    nx = numx3c_xc3c(isorp,index)
    ny = numy3c_xc3c(isorp,index)
    xxmax = x3cmax_xc3c(isorp,index)
    yymax = y3cmax_xc3c(isorp,index)
    if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
       write (*,*) ' What the heck is going on in trescentros!!! error!!! ' 
       write (*,*) ' x = ', x, ' Max of data = ', xxmax
       write (*,*) ' y = ', y, ' Max of data = ', yymax
       stop
    end if
    do iME = 1, index_max3c(in1,in2)
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(0,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(1,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(2,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(3,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, xc3c_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
      bcnalist(4,iME) = Q_L
    end do
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
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(0,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(1,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(2,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(3,iME) = Q_L
      call interpolate_2d (x, y, kforce, nx, ny, hx, hy, den3_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
      bcnalist(4,iME) = Q_L
    end do
  end if
  cost2 = cost**2
  argument = 1.0d0 - cost2
  if (argument .lt. 1.0d-5) argument = 1.0d-5
  sint = sqrt(argument) 
  if (ntheta .ne. 5) then
      write(*,*) ' ntheta must be 5, but it is ',ntheta
      stop       
  end if
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

