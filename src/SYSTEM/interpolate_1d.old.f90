 subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
  use M_constants, only: wp
  use M_system
  use M_constants
  use M_fdata, only: ind2c, numz2c, z2cmax, splineint_2c, nfofx
  implicit none
  integer, intent(in) :: interaction
  integer, intent(in) :: isub
  integer, intent(in) :: in1
  integer, intent(in) :: in2
  integer, intent(in) :: non2c
  integer, intent(in) :: ioption   ! Derivative or not?
  real(wp), intent(in)  :: xin   ! x
  real(wp), intent(out) :: yout  ! F(x)
  real(wp), intent(out) :: dfdx  ! dF/dx
  real(wp), parameter :: tol=1.0d-5
  real(wp), parameter :: e6t=.166666667d0
  real(wp), parameter :: e24t=.04166666667d0
  real(wp), parameter :: e5t=.2d0
  real(wp), parameter :: e2t5=.4d0
  integer i
  integer ileft, imid, iright
  integer iprod, isum
  integer j, jj, jxx
  integer k
  integer nn
  integer nnum
  integer norder
  real(wp) h
  real(wp) xmax
  real(wp), parameter :: xmin = 0.0d0
  real(wp) xxp
  real(wp) f0p1, f0p10, f0p2, f0p3, f0p30, f0p6, f0p8, f1m1, f1m12
  real(wp) f1m16, f1m2, f1m3, f1m4, f1p14, f1p16, f1p2, f1p24, f1p3
  real(wp) f1p4, f1p6, f2m1, f2m3, f2p1, f2p6, f2p7, f3p1, f3p2, ftp
  real(wp) p, pden
  real(wp) prod, prodd 
  real(wp) sum, sumx
  real(wp) xprod, xsumoverj
  real(wp), dimension (5) :: bb
  real(wp), dimension (nfofx) :: pdenom
  real(wp), dimension (nfofx) :: xx
  integer iam
  real(wp), dimension (0:nfofx) :: a, b, c, d
  real(wp), dimension (0:nfofx) :: alpha
  real(wp), dimension (0:nfofx) :: L
  real(wp), dimension (0:nfofx) :: mu
  real(wp), dimension (0:nfofx) :: Z
  real(wp) aaa, bbb, ccc, ddd
  jxx = ind2c(interaction,isub)
  nnum = numz2c(jxx,in1,in2)
  xmax = z2cmax(jxx,in1,in2)
  h = (xmax - xmin)/(nnum - 1)
  if (xin .lt. xmin) then
    write (*,*) ' xin, xmin = ', xin, xmin
    write (*,*) '  error in intrp1d : xin < xmin'
    stop ! Negative value is very bad - should never get there
  else if (xin .gt. xmax) then
    if (abs((xin - xmax)/xmax) .gt. tol) then
     write (*,*) ' xin, xmax = ', xin, xmax
     write (*,*) '  error in intrp1d : xin > xmax'
     write(*,*) 'Interaction = ', interaction
     stop
    end if
    xxp = xmax - xmin
  else if (xin .eq. xmin) then
    yout = splineint_2c(1,non2c,1,jxx,in1,in2)
    if (ioption .eq. 1) dfdx = 0
    return  
  else
    xxp = xin - xmin
  end if
  imid = int(xxp/h) + 1
  if (imid .gt. nnum) imid=nnum ! If we have gone off of the end
  aaa=splineint_2c(1,non2c,imid,jxx,in1,in2)
  bbb=splineint_2c(2,non2c,imid,jxx,in1,in2)
  ccc=splineint_2c(3,non2c,imid,jxx,in1,in2)
  ddd=splineint_2c(4,non2c,imid,jxx,in1,in2)

  xxp=(xxp-(imid-1)*h)/h

  if(ioption .eq. 1) dfdx=(bbb+2.0d0*ccc*xxp+3.0d0*ddd*xxp**2)/h
  yout=aaa+bbb*xxp+ccc*xxp**2+ddd*xxp**3
  return
end

