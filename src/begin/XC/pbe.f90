! $Header: /cvs/begin/XC/pbe.f,v 5.0 2005/01/10 21:07:41 hw93 Exp $
!==========================================================================
! c original version by K Burke
! Perdew-Burke-Ernzerhof GGA
! see: J.P. Perdew, K. Burke, M. Ernzerhof, Phys Rev Lett 77, 3865 (1996).
! this collection contains the older PW91 and Becke exchange as well,
! everything not needed for the PBE GGA is commented out w/ "c--"
!
! Martin Fuchs, FHI der MPG, Berlin, 11-1996
!==========================================================================
!
!--c PBE alpha2.1:
!--c Perdew-Burke-Ernzerhof generalized gradient approximation to the
!--c density functional for exchange-correlation energy of a many-electron
!--c system.
!--c  --------------------------------------------------------------------
!--c |WARNING!  PBE is a simplification of PW91, which yields almost      |
!--c |identical numerical results with simpler formulas from a simpler    |
!--c |derivation.  If you should find significant DIFFERENCES between     |
!--c |PBE and PW91 results, please consult kieron@merlin.phy.tulane.edu   |
!--c |or perdew@mailhost.tcs.tulane.edu.  Thank you.                      |
!--c  --------------------------------------------------------------------
!--c Note: Neglects small grad (zeta) contributions to the correlation
!--c energy.
!--c
!--c Programs implement functional in PBE paper, July 1996 version.
!--c
!--c----------------------------------------------------------------------
!--c Main program testing PBE subroutines for exchange-correlation energy
!--c and potentials, by application to unequal exponential
!--c up- and down-spin densities, showing that the functional derivative
!--c (exchange-correlation potential) correctly generates the energy change
!--c due to a small change in the density.
!--c Kieron Burke, July 2, 1996.
!--C Atomic units are used, so all energies are in hartrees and all
!--c distances in bohrs.
!--c 1 hartree=27.2116eV=627.5kcal/mol; 1bohr=0.529E-10m.
!--c The output should be:
!--c Fup Fdn Zup Zdn             Exc           CHNG1          CHNG
!--c 1.0  .0 1.0  .5  -.311916530229   .000000000000   .0000000000
!--c 1.0  .2 1.0  .5  -.336377065446  -.053880102756  -.0538804290
!--c 1.0  .4 1.0  .5  -.369084886012  -.120463328976  -.1204642921
!--c 1.0  .6 1.0  .5  -.406088525151  -.193370595518  -.1933723422
!--c 1.0  .8 1.0  .5  -.446305936853  -.271252139632  -.2712547575
!--c 1.0 1.0 1.0  .5  -.489150144888  -.353405855349  -.3534094042
!--c 1.0 1.0  .5  .5  -.341059977353  -.316599687356  -.3166037653
!--c 1.0 1.0 1.0 1.0  -.653407740519  -.309758886707  -.3097606837
!--c 1.0 1.0 1.5 1.5  -.962039224827  -.307820467953  -.3078216918
!--c 1.0 1.0 2.0 2.0 -1.269410948459  -.307021487395  -.3070225637
!--c----------------------------------------------------------------------
!--c----------------------------------------------------------------------
!--      IMPLICIT REAL*8 (A-H,O-Z)
!--      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
!--      pi=4.d0*datan(1.d0)
!--      CONF=(3.D0*PI*pi)**THRD
!--      CONRS=(3.D0/(4.D0*PI))**THRD
!--      write(6,*)'Fup Fdn Zup Zdn             Exc'
!--     1,'           CHNG1          CHNG'
!--c----------------------------------------------------------------------
!--c----------------------------------------------------------------------
!--C BEGIN THE LOOP THAT SELECTS A TRIAL DENSITY
!--c spin-densities are of the form
!--c          rho(r)=f*(Z**3/pi)*dexp(-2*Z*r)
!--c delzdn=small change in zdn to test potentials
!--c jdens=counter for which density
!--      DO JDENS = 1,10
!--        FUP=1.D0
!--        FDN=0.2D0*(JDENS-1)
!--        ZUP=1.D0
!--        ZDN=0.5D0
!--        IF(JDENS.GT.6)then
!--       FDN=1.D0
!--          ZUP=0.5D0+0.5D0*(JDENS-7)
!--          ZDN=ZUP
!--     endif
!--        DELZDN=1D-5
!------------------------------------------------------------------
!--C BEGIN THE LOOP THAT INCREMENTS THE DENSITY DIFFERENTIALLY
!--c kdif=1=>density as above
!--c kdif=2=>Zdn changed by DELZDN
!--        DO KDIF=1,2
!--          IF(KDIF.EQ.2)ZDN=ZDN+DELZDN
!------------------------------------------------------------------
!------------------------------------------------------------------
!--C BEGIN THE RADIAL LOOP
!--c sumexc=integrated exchange-correlation energy
!--c chng1=integrated xc energy change, based on vxc
!--c nr=number of points in radial loop
!--c rf=final value of r in integrals
!--c dr=change in r
!--c wt=weight of r in trapezoidal rule
!--c dup=up density
!--c agrup=|grad up|
!--c delgrup=(grad up).(grad |grad up|)
!--c uplap=grad^2 up=Laplacian of up
!--c dn,agrdn,delgrdn,dnlap=corresponding down quantities
!--c d=up+dn
!--c agrad=|grad rho|
!--c delgrad=(grad rho).(grad |grad rho|)
!--          sumexc=0.0D0
!--          CHNG1=0.0D0
!--       nr=10000
!--       rf=20.d0
!--       dr=rf/real(nr)
!--          DO I=1,nr
!--            R=I*dr
!--            WT=4.d0*PI*R*R*dr
!--            DUP=FUP*(ZUP**3/PI)*DEXP(-2.D0*ZUP*R)
!--            DDN=FDN*(ZDN**3/PI)*DEXP(-2.D0*ZDN*R)
!--            ZDNNU=ZDN+DELZDN
!--            DELDDN=FDN*(ZDNNU**3/PI)*DEXP(-2.D0*ZDNNU*R)-DDN
!--         agrup=2.d0*zup*dup
!--         delgrup=8.d0*(zup**3)*dup*dup
!--         uplap=4.d0*zup*dup*(zup-1.d0/r)
!--         agrdn=2.d0*zdn*ddn
!--         delgrdn=8.d0*(zdn**3)*ddn*ddn
!--         dnlap=4.d0*zdn*ddn*(zdn-1.d0/r)
!--            D=DUP+DDN
!--            agrad=2.d0*(ZUP*DUP+ZDN*DDN)
!--         delgrad=4.d0*agrad*(ZUP**2*DUP+ZDN**2*DDN)
!--            call easypbe(dup,agrup,delgrup,uplap,ddn,agrdn,delgrdn,
!--     1           dnlap,agrad,delgrad,1,1,
!--     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
!--     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
!--     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
!--         sumexc=sumexc+d*(expbe+ecpbe)*wt
!--            CHNG1=CHNG1+(vxdnpbe+vcdnpbe)*DELDDN*WT/DELZDN
!--       enddo
!--          IF(KDIF.EQ.1)then
!--         sumEXCO=sumEXC
!--       endif
!--        enddo
!-----------------------------------------------------------------------
!--c----------------------------------------------------------------------
!--C  CHNG: DIRECT XC ENERGY INCREMENT
!--C  IF THE FUNCTIONAL DERIVATIVE IS CORRECT, THEN CHNG1=CHNG
!--        CHNG=(sumEXC-sumEXCO)/DELZDN
!--        PRINT 200,FUP,FDN,ZUP,ZDN,sumEXC,CHNG1,chng
!--      enddo
!--      STOP
!--  200 FORMAT(4f4.1,2f16.12,f14.10)
!--      END
!--c----------------------------------------------------------------------
!--c###############################################################
!--c---------------------------------------------------------------
!--      subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
!--     1           agr,delgr,lcor,lpot,
!--     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
!--     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
!--     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
!--c---------------------------------------------------------------
!--c---------------------------------------------------------------
!-c EASYPBE is a driverforthe PBE subroutines, using simple inputs
!--c K. Burke, May 14, 1996.
!--c inputs: up=up density
!--c    : agrup=|grad up|
!--c    : delgrup=(grad up).(grad |grad up|)
!--c    : uplap=grad^2 up=Laplacian of up
!--c    : dn,agrdn,delgrdn,dnlap=corresponding down quantities
!--c    : agr=|grad rho|
!--c    : delgr=(grad rho).(grad |grad rho|)
!--c    : lcor=flag to do correlation(=0=>don't)
!--c    : lpot=flag to do potential(=0=>don't)
!--c outputs: exlsd=LSD exchange energy density, so that
!--c            ExLSD=int d^3r rho(r) exlsd(r)
!--c     : vxuplsd=up LSD exchange potential
!--c     : vxdnlsd=down LSD exchange potential
!--c        : exclsd=LSD exchange-correlation energy density
!--c     : vxcuplsd=up LSD exchange-correlation potential
!--c     : vxcdnlsd=down LSD exchange-correlation potential
!--c        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
!--c        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
!--c---------------------------------------------------------------
!--c---------------------------------------------------------------
!--c needed constants:
!--c pi32=3 pi**2
!--c alpha=(9pi/4)**thrd
!--      implicit real*8(a-h,o-z)
!--      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
!--      parameter(pi32=29.608813203268075856503472999628d0)
!--      parameter(pi=3.1415926535897932384626433832795d0)
!--      parameter(alpha=1.91915829267751300662482032624669d0)
!--c---------------------------------------------------------------
!--c---------------------------------------------------------------
!--c PBE exchange
!--c use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
!--c do up exchange
!--c fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3)
!--c s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
!--c u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
!--c v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
!--      rho2=2.d0*up
!--      if(rho2.gt.1.d-15)then
!--        fk=(pi32*rho2)**thrd
!--        s=2.d0*agrup/(2.d0*fk*rho2)
!--        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
!--        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
!--        call exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
!--        call exchpw91(rho2,s,u,v,exuppw91,vxuppw91)
!--        call exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
!--      else
!--     exuplsd=0.d0
!--     vxuplsd=0.d0
!--     exuppw91=0.d0
!--     vxuppw91=0.d0
!--     exuppbe=0.d0
!--     vxuppbe=0.d0
!--      endif
!--c repeat for down
!--      rho2=2.d0*dn
!--      if(rho2.gt.1.d-15)then
!--        fk=(pi32*rho2)**thrd
!--        s=2.d0*agrdn/(2.d0*fk*rho2)
!--        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
!--        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
!--        call exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
!--        call exchpw91(rho2,s,u,v,exdnpw91,vxdnpw91)
!--        call exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
!--      else
!--     exdnlsd=0.d0
!--     vxdnlsd=0.d0
!--     exdnpw91=0.d0
!--     vxdnpw91=0.d0
!--     exdnpbe=0.d0
!--     vxdnpbe=0.d0
!--      endif
!--10    continue
!--c construct total density and contribution to ex
!--      rho=up+dn
!--      if(rho .gt. 1.d-15) then
!--        exlsd=(exuplsd*up+exdnlsd*dn)/rho
!--        expw91=(exuppw91*up+exdnpw91*dn)/rho
!--        expbe=(exuppbe*up+exdnpbe*dn)/rho
!--      else
!--        exlsd=0.d0
!--        expw91=0.d0
!--        expbe=0.d0
!--      endif
!--      if(lcor.eq.0)return
!--c---------------------------------------------------------------
!--c---------------------------------------------------------------
!--c Now do correlation
!--c zet=(up-dn)/rho
!--c g=phi(zeta)
!--c rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
!--c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
!--c twoksg=2*Ks*phi
!--c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
!--c uu=delgrad/(rho^2*twoksg^3)
!--c rholap=Laplacian
!--c vv=Laplacian/(rho*twoksg^2)
!--c ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
!--c ec=lsd correlation energy
!--c vcup=lsd up correlation potential
!--c vcdn=lsd down correlation potential
!--c h=gradient correction to correlation energy
!--c dvcup=gradient correction to up correlation potential
!--c dvcdn=gradient correction to down correlation potential
!--      if(rho.lt.1.d-18)return
!--      zet=(up-dn)/rho
!--      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
!--      fk=(pi32*rho)**thrd
!--      rs=alpha/fk
!--      sk=sqrt(4.d0*fk/pi)
!--      twoksg=2.d0*sk*g
!--      t=agr/(twoksg*rho)
!--      uu=delgr/(rho*rho*twoksg**3)
!--      rholap=uplap+dnlap
!--      vv=rholap/(rho*twoksg**2)
!--      ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
!--      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn,
!--     1                  H,DVCUP,DVCDN)
!--      eclsd=ec
!--      ecpbe=ec+h
!--      vcuplsd=vcup
!--      vcdnlsd=vcdn
!--      vcuppbe=vcup+dvcup
!--      vcdnpbe=vcdn+dvcdn
!--      call CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!--      call CORPW91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,DVCUP,DVCDN)
!--      ecpw91=ec+h
!--      vcuppw91=vcup+dvcup
!--      vcdnpw91=vcdn+dvcdn
!--      return
!--      end
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke's modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  input lgga:  (=0=>don't put in gradient corrections, just LDA)
!  input lpot:  (=0=>don't get potential and don't need U and V)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      use precision
      implicit none
      real(kind=long), intent(in):: rho, s, u, v, lgga, lpot
      real(kind=long), intent(out):: ex, vx
      real(kind=long), parameter :: thrd = 0.333333333333d0
      real(kind=long), parameter :: thrd4 = 1.333333333333d0
      real(kind=long), parameter :: ax = -0.738558766382022405884230032680836d0
      real(kind=long), parameter :: um = 0.2195149727645171d0
      real(kind=long), parameter :: uk = 0.8040d0
      real(kind=long) :: ul, exunif, s2, p0, fxpbe, fs, fss
!     parameter(pi=3.14159265358979323846264338327950d0)
      ul = um/uk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
        ex=exunif
        vx=ex*thrd4
        return
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24)
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,H,DVCUP,DVCDN)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
!       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!       :  UU,VV,WW, only needed for PBE potential
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      use precision
      implicit none
      real(kind=long), intent(in) :: rs, zet, t, uu, vv, ww, lgga, lpot
      real(kind=long), intent(out) :: ec, vcup, vcdn, h, dvcup, dvcdn
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
      real(kind=long), parameter :: thrd = 0.333333333333d0
      real(kind=long), parameter :: thrdm = -0.333333333333d0
      real(kind=long), parameter :: thrd2 = 0.666666666667d0
      real(kind=long), parameter :: sixthm = -0.16666666667d0
      real(kind=long), parameter :: thrd4 = 1.33333333333d0
      real(kind=long), parameter :: GAM = 0.5198420997897463295344212145565d0
      real(kind=long), parameter :: fzz = 1.70992093416d0
      real(kind=long), parameter :: gamma = 0.03109069086965489503494086371273d0
      real(kind=long), parameter :: bet = 0.06672455060314922d0
      real(kind=long), parameter :: delt = 2.14612633997d0
      real(kind=long), parameter :: eta = 1.0d-12

      real(kind=long) :: EU
      real(kind=long) :: EURS
      real(kind=long) :: EP
      real(kind=long) :: EPRS
      real(kind=long) :: ALFM
      real(kind=long) :: ALFRSM
      real(kind=long) :: rtrs
      real(kind=long) :: ALFC
      real(kind=long) :: Z4
      real(kind=long) :: F
      real(kind=long) :: ECRS
      real(kind=long) :: FZ
      real(kind=long) :: ECZET
      real(kind=long) :: COMM
      real(kind=long) :: G
      real(kind=long) :: G3
      real(kind=long) :: PON
      real(kind=long) :: B
      real(kind=long) :: B2
      real(kind=long) :: T2
      real(kind=long) :: T4
      real(kind=long) :: RS2
      real(kind=long) :: RS3
      real(kind=long) :: Q4
      real(kind=long) :: Q5
      real(kind=long) :: G4
      real(kind=long) :: T6
      real(kind=long) :: RSTHRD
      real(kind=long) :: GZ
      real(kind=long) :: FAC
      real(kind=long) :: BG
      real(kind=long) :: BEC
      real(kind=long) :: Q8
      real(kind=long) :: Q9
      real(kind=long) :: hB
      real(kind=long) :: hRS
      real(kind=long) :: FACT0
      real(kind=long) :: FACT1
      real(kind=long) :: hBT
      real(kind=long) :: hRST
      real(kind=long) :: hZ
      real(kind=long) :: hT
      real(kind=long) :: hZT
      real(kind=long) :: FACT2
      real(kind=long) :: FACT3
      real(kind=long) :: hTT
      real(kind=long) :: PREF
      real(kind=long) :: FACT5

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=sqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(EXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*LOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm-((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      use precision
      implicit none
      real(kind=long), intent(in) :: a, a1, b1, b2, b3, b4, rtrs
      real(kind=long), intent(out) :: gg, ggrs
      real(kind=long) :: q0, q1, q2, q3
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
!--c----------------------------------------------------------------------
!--c######################################################################
!--c----------------------------------------------------------------------
!--      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
!--C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!--C  INPUT D : DENSITY
!--C  INPUT S:  ABS(GRAD D)/(2*KF*D)
!--C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
!--C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
!--C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!--      IMPLICIT REAL*8 (A-H,O-Z)
!--      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
!--      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
!--      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
!--c for Becke exchange, set a3=b1=0
!--      FAC = AX*D**THRD
!--      S2 = S*S
!--      S3 = S2*S
!--      S4 = S2*S2
!--      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
!--      P1 = DLOG(A*S+1.D0/P0)
!--      P2 = DEXP(-A4*S2)
!--      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
!--      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
!--      F = P3*P4
!--      EX = FAC*F
!--C  LOCAL EXCHANGE OPTION
!--C     EX = FAC
!--C  ENERGY DONE. NOW THE POTENTIAL:
!--      P5 = B1*S2-(A2-A3*P2)
!--      P6 = A1*S*(P1+A*S*P0)
!--      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
!--      FS = P3*(P3*P5*P6+P7)
!--      P8 = 2.D0*S*(B1-A3*A4*P2)
!--      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
!--      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
!--      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
!--      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
!--      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
!--C  LOCAL EXCHANGE OPTION:
!--C     VX = FAC*THRD4
!--      RETURN
!--      END
!--c----------------------------------------------------------------------
!--c######################################################################
!--c----------------------------------------------------------------------
!--      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!--C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!--C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!--C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!--C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!--C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!--      IMPLICIT REAL*8 (A-H,O-Z)
!--      parameter(gam=0.5198421D0,fzz=1.709921D0)
!--      parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
!--      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
!--      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
!--     1    0.49294D0,1.00D0,RS,EU,EURS)
!--      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
!--     1    0.62517D0,1.00D0,RS,EP,EPRS)
!--      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
!--     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
!--C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
!--      ALFC = -ALFM
!--      Z4 = ZET**4
!--      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!--C  ENERGY DONE. NOW THE POTENTIAL:
!--      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
!--      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
!--      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
!--     1        -(1.D0-Z4)*ALFM/FZZ)
!--      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
!--      VCUP = COMM + ECZET
!--      VCDN = COMM - ECZET
!--      RETURN
!--      END
!--c----------------------------------------------------------------------
!--c######################################################################
!--c----------------------------------------------------------------------
!--      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
!--C  CALLED BY SUBROUTINE CORLSD
!--      IMPLICIT REAL*8 (A-H,O-Z)
!--      P1 = P + 1.D0
!--      Q0 = -2.D0*A*(1.D0+A1*RS)
!--      RS12 = DSQRT(RS)
!--      RS32 = RS12**3
!--      RSP = RS**P
!--      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
!--      Q2 = DLOG(1.D0+1.D0/Q1)
!--      GG = Q0*Q2
!--      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
!--      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
!--      RETURN
!--      END
!--c----------------------------------------------------------------------
!--c######################################################################
!--c----------------------------------------------------------------------
!--      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
!--     1                   DVCUP,DVCDN)
!--C  pw91 CORRELATION, modified by K. Burke to put all arguments
!--c  as variables in calling statement, rather than in common block
!--c  May, 1996.
!--C  INPUT RS: SEITZ RADIUS
!--C  INPUT ZET: RELATIVE SPIN POLARIZATION
!--C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
!--C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
!--C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
!--C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
!--C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!--C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
!--      IMPLICIT REAL*8 (A-H,O-Z)
!--      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
!--      parameter(alf=0.09D0)
!--      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
!--      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
!--      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
!--      BET = XNU*CC0
!--      DELT = 2.D0*ALF/BET
!--      G3 = G**3
!--      G4 = G3*G
!--      PON = -DELT*EC/(G3*BET)
!--      B = DELT/(DEXP(PON)-1.D0)
!--      B2 = B*B
!--      T2 = T*T
!--      T4 = T2*T2
!--      T6 = T4*T2
!--      RS2 = RS*RS
!--      RS3 = RS2*RS
!--      Q4 = 1.D0+B*T2
!--      Q5 = 1.D0+B*T2+B2*T4
!--      Q6 = C1+C2*RS+C3*RS2
!--      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
!--      CC = -CX + Q6/Q7
!--      R0 = 0.663436444d0*rs
!--      R1 = A4*R0*G4
!--      COEFF = CC-CC0-3.D0*CX/7.D0
!--      R2 = XNU*COEFF*G3
!--      R3 = DEXP(-R1*T2)
!--      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
!--      H1 = R3*R2*T2
!--      H = H0+H1
!--C  LOCAL CORRELATION OPTION:
!--C     H = 0.0D0
!--C  ENERGY DONE. NOW THE POTENTIAL:
!--      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
!--      RSTHRD = RS/3.D0
!--      R4 = RSTHRD*CCRS/COEFF
!--      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
!--      FAC = DELT/B+1.D0
!--      BG = -3.D0*B2*EC*FAC/(BET*G4)
!--      BEC = B2*FAC/(BET*G3)
!--      Q8 = Q5*Q5+DELT*Q4*Q5*T2
!--      Q9 = 1.D0+2.D0*B*T2
!--      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
!--      H0RS = -RSTHRD*H0B*BEC*ECRS
!--      FACT0 = 2.D0*DELT-6.D0*B
!--      FACT1 = Q5*Q9+Q4*Q9*Q9
!--      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
!--      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
!--      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
!--      H0T = 2.*BET*G3*Q9/Q8
!--      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
!--      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
!--      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
!--      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
!--      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
!--      FACT4 = 2.D0-R1*T2
!--      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
!--      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
!--      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
!--      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
!--      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
!--      HRS = H0RS+H1RS
!--      HRST = H0RST+H1RST
!--      HT = H0T+H1T
!--      HTT = H0TT+H1TT
!--      HZ = H0Z+H1Z
!--      HZT = H0ZT+H1ZT
!--      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
!--      PREF = HZ-GZ*T2*HT/G
!--      FACT5 = GZ*(2.D0*HT+T*HTT)/G
!--      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
!--      DVCUP = COMM + PREF
!--      DVCDN = COMM - PREF
!--C  LOCAL CORRELATION OPTION:
!--C     DVCUP = 0.0D0
!--C     DVCDN = 0.0D0
!--      RETURN
!--      END
!--c----------------------------------------------------------------------
