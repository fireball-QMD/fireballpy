 subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
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
  real*8, intent(in)  :: xin   ! x
  real*8, intent(out) :: yout  ! F(x)
  real*8, intent(out) :: dfdx  ! dF/dx
  real*8, parameter :: tol=1.0d-5
  real*8, parameter :: e6t=.166666667d0
  real*8, parameter :: e24t=.04166666667d0
  real*8, parameter :: e5t=.2d0
  real*8, parameter :: e2t5=.4d0
  integer i
  integer ileft, imid, iright
  integer iprod, isum
  integer j, jj, jxx
  integer k
  integer nn
  integer nnum
  integer norder
  real*8 h
  real*8 xmax
  real*8, parameter :: xmin = 0.0d0
  real*8 xxp
  real*8 f0p1, f0p10, f0p2, f0p3, f0p30, f0p6, f0p8, f1m1, f1m12
  real*8 f1m16, f1m2, f1m3, f1m4, f1p14, f1p16, f1p2, f1p24, f1p3
  real*8 f1p4, f1p6, f2m1, f2m3, f2p1, f2p6, f2p7, f3p1, f3p2, ftp
  real*8 p, pden
  real*8 prod, prodd 
  real*8 sum, sumx
  real*8 xprod, xsumoverj
  real*8, dimension (5) :: bb
  real*8, dimension (nfofx) :: pdenom
  real*8, dimension (nfofx) :: xx
  integer iam
  real*8, dimension (0:nfofx) :: a, b, c, d
  real*8, dimension (0:nfofx) :: alpha
  real*8, dimension (0:nfofx) :: L
  real*8, dimension (0:nfofx) :: mu
  real*8, dimension (0:nfofx) :: Z
  real*8 aaa, bbb, ccc, ddd
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

