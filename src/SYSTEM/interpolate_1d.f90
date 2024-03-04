 subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
  use M_system
  use M_constants
  use M_fdata, only: ind2c, numz2c, z2cmax, xintegral_2c, nfofx
  implicit none
  integer, intent(in) :: interaction
  integer, intent(in) :: isub
  integer, intent(in) :: in1
  integer, intent(in) :: in2
  integer, intent(in) :: non2c
  integer, intent(in) :: ioption   ! Derivative or not?
  real, intent(in)  :: xin   ! x
  real, intent(out) :: yout  ! F(x)
  real, intent(out) :: dfdx  ! dF/dx
  real, parameter :: tol=1.0d-5
  real, parameter :: e6t=.166666667d0
  real, parameter :: e24t=.04166666667d0
  real, parameter :: e5t=.2d0
  real, parameter :: e2t5=.4d0
  integer i
  integer ileft, imid, iright
  integer iprod, isum
  integer j, jj, jxx
  integer k
  integer nn
  integer nnum
  integer norder
  real h
  real xmax
  real, parameter :: xmin = 0.0d0
  real xxp
  real f0p1, f0p10, f0p2, f0p3, f0p30, f0p6, f0p8, f1m1, f1m12
  real f1m16, f1m2, f1m3, f1m4, f1p14, f1p16, f1p2, f1p24, f1p3
  real f1p4, f1p6, f2m1, f2m3, f2p1, f2p6, f2p7, f3p1, f3p2, ftp
  real p, pden
  real prod, prodd 
  real sum, sumx
  real xprod, xsumoverj
  real, dimension (5) :: bb
  real, dimension (nfofx) :: pdenom
  real, dimension (nfofx) :: xx
  integer iam
  real, dimension (0:nfofx) :: a, b, c, d
  real, dimension (0:nfofx) :: alpha
  real, dimension (0:nfofx) :: L
  real, dimension (0:nfofx) :: mu
  real, dimension (0:nfofx) :: Z
  real aaa, bbb, ccc, ddd
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
    yout = xintegral_2c(non2c,1,jxx,in1,in2)
    if (ioption .eq. 1) dfdx = 0
    return  
  else
    xxp = xin - xmin
  end if
  imid = int(xxp/h) + 1
  if (imid .gt. nnum) imid=nnum ! If we have gone off of the end
  norder=abs(norder1)
  if(norder .eq. 0 .or. norder+1 .gt. nnum)then
    write(*,*) ' norder1 is chosen wrong in interpolate_1d'
    write(*,*) ' norder=',norder1,' and must be at least +/- 1'
    stop
  end if
  if(norder1 .eq. 3) go to 93
  if(norder1 .eq. 5) go to 96
  if(mod(norder+1,2).eq.0)then
     ileft=imid-((norder-1)/2)
     iright=imid+((norder+1)/2)
  else
     ileft=imid-(norder/2)
     iright=imid+(norder/2)
  endif
  if(ileft.lt.1)then
     ileft=1
     iright=norder+1
  elseif(iright.gt.nnum)then
     ileft=nnum-norder
     iright=nnum
  endif
  if(norder1 .lt. 0) goto 111
  do i=ileft,iright
   xx(i)=(i-1)*h
  end do
  sum=0.0e0
  do isum=ileft,iright
    prod=1.0e0
    pdenom(isum)=1.0e0
    do iprod=ileft,iright
      if(iprod.ne.isum)then
        pden=1.0e0/(xx(isum)-xx(iprod))
        pdenom(isum)=pdenom(isum)*pden
        prod=prod*(xxp-xx(iprod))*pden
      end if
    end do
    sum=sum+xintegral_2c(non2c,isum,jxx,in1,in2)*prod
  end do
  yout=sum
  dfdx=0.0e0
  if(ioption.ne.1) return
  sumx=0.0e0
  do 1208 i=ileft,iright
   xsumoverj=0.0e0
   do 1308 j=ileft,iright
    if(j.eq.i)goto 1308
    xprod=1.0e0
    do 1408 k=ileft,iright
     if(k.eq.j)goto 1408
     if(k.eq.i)goto 1408
     xprod=xprod*(xxp-xx(k))
1408      continue
    xsumoverj=xsumoverj+xprod
1308     continue
   sumx=sumx+xintegral_2c(non2c,i,jxx,in1,in2)*pdenom(i)*xsumoverj
1208    continue
  dfdx=sumx
  return
  93    continue
  if(imid .lt. 2) imid=2
  nn=nnum-2
  if(imid .gt. nn) imid=nn
  ftp=xintegral_2c(non2c,imid-1,jxx,in1,in2)
  f1m1=ftp
  f1m2=ftp*2
  f1m3=ftp*3
  ftp=xintegral_2c(non2c,imid,jxx,in1,in2)
  f0p1=ftp
  f0p3=ftp*3
  f0p6=ftp*6
  ftp=xintegral_2c(non2c,imid+1,jxx,in1,in2)
  f1p3=ftp*3
  f1p6=ftp*6
  f2p1=xintegral_2c(non2c,imid+2,jxx,in1,in2)
  bb(3)=-f1m1+f0p3-f1p3+f2p1
  bb(2)=f1m3-f0p6+f1p3
  bb(1)=-f1m2-f0p3+f1p6-f2p1
  p=(xxp-h*(imid-1))/h
  prod=bb(3)*p
  prodd=3*prod
  do jj=2,1,-1
   ftp=bb(jj)
   prod=(prod+ftp)*p
   prodd=prodd+jj*ftp
   if(jj .ne. 1) prodd=prodd*p
  end do
  yout=e6t*prod+f0p1
  dfdx=e6t*prodd/h
  return
  96    continue
  if(imid .lt. 3) imid=3
  nn=nnum-3
  if(imid .gt. nn) imid=nn
  ftp=xintegral_2c(non2c,imid-2,jxx,in1,in2)
  f2m1=ftp
  f2m3=ftp*3
  ftp=xintegral_2c(non2c,imid-1,jxx,in1,in2)
  f1m1=ftp
  f1m4=4*ftp
  f1m16=16*ftp
  f1m12=12*ftp
  ftp=xintegral_2c(non2c,imid,jxx,in1,in2)
  f0p1=ftp
  f0p2=2*ftp
  f0p6=6*ftp
  f0p8=8*ftp
  f0p10=10*ftp
  f0p30=30*ftp
  ftp=xintegral_2c(non2c,imid+1,jxx,in1,in2)
  f1p2=2*ftp
  f1p4=4*ftp
  f1p14=14*ftp
  f1p16=16*ftp
  f1p24=24*ftp
  ftp=xintegral_2c(non2c,imid+2,jxx,in1,in2)
  f2p1=ftp
  f2p6=6*ftp
  f2p7=7*ftp
  ftp=xintegral_2c(non2c,imid+3,jxx,in1,in2)
  f3p1=ftp
  f3p2=ftp*2
  bb(5)= f1m1-f0p2+f1p2-f2p1+e5t*(f3p1-f2m1)
  bb(4)= f2m1-f1m4+f0p6-f1p4+f2p1
  bb(3)=-f2m1-f1m1+f0p10-f1p14+f2p7-f3p1
  bb(2)=-f2m1+f1m16-f0p30+f1p16-f2p1
  bb(1)=-f1m12-f0p8+f1p24-f2p6+e2t5*(f3p2+f2m3)
  p=(xxp-h*(imid-1))/h
  prod=bb(5)*p
  prodd=prod*5
  do jj=4,1,-1
    ftp=bb(jj)
    prod=(prod+ftp)*p
    prodd=prodd+jj*ftp
    if(jj .ne. 1) prodd=prodd*p
  end do
  yout=e24t*prod+f0p1
  dfdx=e24t*prodd/h
  return
 111    continue
  do i=0,norder
    a(i)=xintegral_2c(non2c,i+ileft,jxx,in1,in2)
  end do
  do i=1,norder-1
    alpha(i)=3.0d0*(a(i+1)-2*a(i)+a(i-1))/h
  end do
  L(0)=1
  mu(0)=0
  Z(0)=0
  c(0)=0
  do i=1,norder-1
    L(i)=(4.0d0-mu(i-1))*h
    mu(i)=h/L(i)
    Z(i)=(alpha(i)-h*Z(i-1))/L(i)
  end do
  L(norder)=1
  mu(norder)=0
  Z(norder)=0
  c(norder)=0
  if(imid .eq. nnum)imid=nnum-1
  iam=imid-ileft
  do j=norder-1,iam,-1
    c(j)=z(j)-mu(j)*c(j+1)
    b(j)=(a(j+1)-a(j))/h-h*(c(j+1)+2.0d0*c(j))/3.0d0
    d(j)=(c(j+1)-c(j))/(3.0d0*h)
  end do
  xxp=xxp-(imid-1)*h
  yout=a(iam)+b(iam)*xxp+c(iam)*xxp**2+d(iam)*xxp**3
  if (ioption .eq. 1) dfdx=b(iam)+2.0d0*c(iam)*xxp+3.0d0*d(iam)*xxp**2
  return
end

