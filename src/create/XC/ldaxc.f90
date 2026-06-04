!**********************************************************************
!  LDA exchange
!  Hartree a.u.
!***********************************************************************
      subroutine xlda(d,vx,ex)
!
      use precision, only: wp
      implicit none

      real(kind=wp) AX,D,EPS,EX,THD,THD4,VX

      data ax,thd,thd4,eps/-.738558766382022406d0 &
     &,                     .333333333333333333d0 &
     &,                     .133333333333333333d1 &
     &,                    1.d-100/
!
      if(d .le. eps) then
        vx=0.d0
        ex=0.d0
      else
        ex=ax*d**thd
        vx=thd4*ex
      endif
!
      return
      end
!
!***********************************************************************
! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).
! Hartree a.u.
! see W.E. Pickett, Comp.Phys.Rep. 9, 115 (1989), pg. 187 : w2 = 7.79
!***********************************************************************
!
      subroutine wigner(rh,ex,fx,exc,fxc)
      use precision, only: wp
      implicit none

      real(kind=wp) AX,EPS,EX,EXC,FX,FXC,PIF,RH,RS,THRD,W1,W2,X,Y

      data thrd,w1,w2,ax,pif,eps/.333333333333333333d0 &
     &,                         -.44d0 &
     &,                          .779d1 &
     &,                         -.738558766382022447d0 &
     &,                          .620350490899400087d0 &
     &,                         1.d-100/

!
      if(rh .lt. eps) then
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
        x  = rh**thrd
        ex = ax * x
        fx = 4.d0*thrd*ex
        rs = pif/x
        y  = 1.d0/(rs + w2)
        exc= ex + w1*y
        fxc= fx + w1*(4.d0*thrd*rs + w2)*y*y
      endif
!
      return
      end

!
!c**********************************************************************
! LDA - scaled Wigner exchange-correlation
! from [Q. Zhao, R.G. Parr, PRA 46, R5320 (1992)] eqs. (5) and (7)
!
! Input
! rh        density
!
! Output
! ex        exchange energy per electron
! fx           ""    potential
!           ... based on continuum-LDA (electron-gas like)
!
! exc       exchange-correlation energy per electron
! fxc          ""        ""      potential
!
! Hartree a.u.
! Martin Fuchs, FHI der MPG, Berlin, 01-1993
!c**********************************************************************
      subroutine wigscaled(rh,ex,fx,exc,fxc)
      use precision, only: wp
      implicit none

      real(kind=wp) AA,AX,BB,EX,EXC,FX,FXC,RH,THRD,X,XYLN,Y,Z

      data thrd,aa,bb,ax    /.333333333333333333d0 &
     &,                      .93222d0 &
     &,                      .947362d-2 &
     &,                      .738558766382022406d0/
!
      if(rh .le. 0.d0) then
        rh = 0.d0
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
        x = bb*rh**thrd
        y = x/(x + 1.d0)
        z = -aa*x/bb
        xyln = x*log(y)
!
        exc = z * (1.d0 + xyln)
        fxc = thrd*z *(4.d0 + 5.d0*xyln + y)
!
! electron-gas based LDA exchange
!
        ex = ax*z/aa
        fx = 4.d0*thrd*ex
      endif
!
      return
      end
!