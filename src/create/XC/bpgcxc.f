!
!  Becke exchange for a spin-unpolarized electronic system
!
!  Gradient-corrected exchange energy based on
!     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
!  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
!     and the functional derivative formula are given by
!     [J.P. Perdew , PRB 33, 8800, 1986].
!     [J.P. Perdew , PRB 34, 7406, 1986].
!  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
!
!  Hartree a.u.
!
!  input
!  d            density
!  s            abs(grad d)/(2kf*d)
!  u            (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!           >>  grad(abs(grad d) has mixed derivatives ! <<
!  v            (laplacian d)/(d*(2*kf)**2)
!
!  output
!  ex           exchange energy per electron
!  vx           exchange potential
!
! Martin Fuchs, FHI der MPG, Berlin, 02-1993
!

      subroutine xbecke(d,s,u,v,ex,vx)

      implicit none

      real*8 d,s,u,v,ex,vx,c,b,bb,ax,thrd,thrd4
      real*8 dd1,ddi,f,fac,fs,fss,g,g1,x,y0,y1,y2

      data c / .779555417944150792d1/
      data b / .42d-2/
      data bb/-.451357747124625192d-2/
      data ax/-.738558766382022406d0/
      data thrd,thrd4/.333333333333333333d0,.1333333333333333333d1/

! exchange enhancement factor f
      x  = c*s
      y1 = 1.d0/sqrt(1.d0+x*x)
      y0 = log(x+1.d0/y1)
      y2 = -x*y1*y1*y1
      ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
      dd1= 6.d0*b*(y0+x*y1)
      g  = 1.d0 - 0.5d0*x*dd1*ddi
      fs = -2.d0*bb*c*c*ddi
      g1 = -3.d0*b*(y0+x*(3.d0*y1+x*y2-dd1*dd1*ddi/(6.d0*b)))
      fss= fs*c*(g1 - g*dd1)*ddi
      fs = fs*g
      f  = 1.d0 - bb*x*x*ddi

! LDA only
      fac= ax*d**thrd

! energy
      ex = fac*f

! potential
      vx = fac*(thrd4*f-(u-thrd4*s*s*s)*fss-v*fs)

      return
      end

!***********************************************************************

!  gradient-correction to correlation energy from
!
!  [J.P.Perdew, PRB 33, 8822 (1986) and PRB 34, 7406 (1986)]
!
!  to correlation part of the Becke-Perdew gradient-corrected
!  xc-functional
!
!  input
!  d1,d2 : up/down spindensity
!  dp12..: grad(d1)*grad(d2) {* == vector product}, MUST NOT == 0
!  uu    : (grad d)*grad(abs(grad d)) , d = d1 + d2
!  vv    : laplacian d
!
!  output
!  dec  : correction to correlation energy per electron
!  dvcup :         - "" -            potential maj. spin
!  dvcdn :         - "" -            potential min. spin
!
!  Hartree a.u.
!
!  Martin Fuchs, FHI der MPG, Berlin, 07-1993
!
!
      subroutine corga86(d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn)
!
      implicit none

      real*8 d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn
      real*8 a1,a2,a3,a4,a5,a6,a7
      real*8 t13,t23,t43,t53,t76,crt2
      real*8 c,c1,c2,cp,d,d43,dd,dm13,dpnorm,dpnorm2,fi
      real*8 uuu,vvv,www,zzz,zet,zz1,zz2
!
      data a1,a2,a3,a4,a5,a6,a7/ 2.568d-3
     &,                           1.443307452d-2
     &,                           2.843543831d-6
     &,                           5.411317331d0
     &,                           1.816419933d-1
     &,                           1.763993811d-2
     &,                           8.12908d-4/
      data t13,t23,t43,t53,t76/.33333333333333333d0
     &,                        .66666666666666667d0
     &,                        .13333333333333333d1
     &,                        .16666666666666667d1
     &,                        .11666666666666667d1/
      data crt2/.1587401052d1/
!
      d    = d1+d2
      dm13 = 1.d0/d**t13
      d43  = d**t43
!
! gradient expansion coefficient
      c1 = a1 + dm13*(a2+dm13*a3)
      c2 = 1.d0 + dm13*(a4+dm13*(a5+dm13*a6))
      c  = 1.667d-3 + c1/c2
!
      dpnorm = sqrt(dp11+dp22+2.d0*dp12)
      dpnorm2= dpnorm*dpnorm
!
      fi = a7*dpnorm/(c*d**t76)
!
! spin interpolation
      zet= (d1-d2)/d
      dd = sqrt(.5d0*((1.d0+zet)**t53 + (1.d0-zet)**t53))
!
! dC(n)/dn
      cp = -t13/(d43*c2*c2)
     &    *(c2*(a2+2.d0*a3*dm13)-c1*(a4+dm13*(2.d0*a5+dm13*3.d0*a6)))
!
! spin-independent term
      www=( (fi-1.d0)*(cp/c-t43/d)+(fi-2.d0)*fi*(cp/c+t76/d) )*dpnorm2
      uuu=(3.d0-fi)*fi*uu/dpnorm
      vvv=(fi-2.d0)*vv
!
! spin dependent term
      zzz=crt2*.5d0*t53/(dd*d43)**2*(d1**t23-d2**t23)
      zz1= zzz*((1.d0-fi)*d2*dpnorm2-(2.d0-fi)*d*(dp12+dp22))
      zz2=-zzz*((1.d0-fi)*d1*dpnorm2-(2.d0-fi)*d*(dp12+dp11))
!
      dec=c/(dd*exp(fi)*d43)
!
      dvcup=dec*(www+uuu+vvv+zz1)
      dvcdn=dec*(www+uuu+vvv+zz2)
      dec=dec*dpnorm2/d
!
      return
      end
