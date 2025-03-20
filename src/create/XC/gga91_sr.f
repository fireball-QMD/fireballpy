!
!  Perdew - Wang GGA91 exchange-correlation functional
!
!  J.P. Perdew et.al., Phys.Rev.B 46, 6671 (1992)
!
!  w/out gradients it's the new parametrization of the Ceperley-Alder
!  xc-energy data from
!
!  J.P. Perdew et.al., Phys.Rev.B 45, 13244 (1992)
!
!  pure spin singularity numerically removed
!
! original version by J P Perdew
! modified -- M. Fuchs (cmf), FHI Berlin, 07-1992
!
!
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  input u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!  input v: (laplacian d)/(d*(2*kf)**2)
!  output:  exchange energy per electron (ex) and potential (vx)
!
      subroutine exch(d,s,u,v,ex,vx)
!
      implicit none
      real*8 d,s,u,v,ex,vx,a1,a2,a3,a4,ax,a,b1,thrd,thrd4
      real*8 f,fac,fs,fss,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
      real*8 s2,s3,s4

      data a1,a2,a3,a4/0.19645d0,0.27430d0,0.15084d0,100.d0/
      data ax,a,b1/-0.7385588d0,7.7956d0,0.004d0/
      data thrd,thrd4/0.333333333333d0,1.33333333333d0/
      fac = ax*d**thrd
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      p0 = 1.d0/sqrt(1.d0+a*a*s2)
      p1 = log(a*s+1.d0/p0)
      p2 = dexp(-a4*s2)
      p3 = 1.d0/(1.d0+a1*s*p1+b1*s4)
      p4 = 1.d0+a1*s*p1+(a2-a3*p2)*s2
      f = p3*p4
      ex = fac*f
!  local exchange option
!     ex = fac
!  energy done. now the potential:
      p5 = b1*s2-(a2-a3*p2)
      p6 = a1*s*(p1+a*s*p0)
      p7 = 2.d0*(a2-a3*p2)+2.d0*a3*a4*s2*p2-4.d0*b1*s2*f
      fs = p3*(p3*p5*p6+p7)
      p8 = 2.d0*s*(b1-a3*a4*p2)
      p9 = a1*p1+a*a1*s*p0*(3.d0-a*a*s2*p0*p0)
      p10 = 4.d0*a3*a4*s*p2*(2.d0-a4*s2)-8.d0*b1*s*f-4.d0*b1*s3*fs
      p11 = -p3*p3*(a1*p1+a*a1*s*p0+4.d0*b1*s3)
      fss = p3*p3*(p5*p9+p6*p8)+2.d0*p3*p5*p6*p11+p3*p10+p7*p11
      vx = fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)
!
!  local exchange option:
!     vx = fac*thrd4
      return
      end

!******************************************************************

      subroutine corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
!  uniform-gas correlation of perdew and wang 1991
!  input: seitz radius (rs), relative spin polarization (zet)
!  output: correlation energy per electron (ec), up- and down-spin
!  potentials (vcup,vcdn), derivatives of ec wrt rs (ecrs) & zet (eczet)
!  output: correlation contribution (alfc) to the spin stiffness
      implicit none
      real*8 rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc
      real*8 gam,fzz,thrd,thrd4
      real*8 ALFM,ALFRSM,COMM,EP,EPRS,EU,EURS,F,FZ,Z4
      data gam,fzz/0.5198421d0,1.709921d0/
      data thrd,thrd4/0.333333333333d0,1.333333333333d0/
      f = ((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
      call gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,
     1    0.49294d0,1.00d0,rs,eu,eurs)
      call gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,
     1    0.62517d0,1.00d0,rs,ep,eprs)
      call gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,
     1    0.49671d0,1.00d0,rs,alfm,alfrsm)
!  alfm is minus the spin stiffness alfc
      alfc = -alfm
      z4 = zet**4
      ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
!  energy done. now the potential:
      ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
      fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
      eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu
     1        -(1.d0-z4)*alfm/fzz)
      comm = ec -rs*ecrs/3.d0-zet*eczet
      vcup = comm + eczet
      vcdn = comm - eczet
      return
      end

!********************************************************************

      subroutine gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
!
      implicit none
      real*8 a,a1,b1,b2,b3,b4,p,rs,gg,ggrs,p1,q0,rs12,rs32,rsp
      real*8 q1,q2,q3
      p1 = p + 1.d0
      q0 = -2.d0*a*(1.d0+a1*rs)
      rs12 = sqrt(rs)
      rs32 = rs12**3
      rsp = rs**p
      q1 = 2.d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2 = log(1.d0+1.d0/q1)
      gg = q0*q2
      q3 = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*p1*rsp)
      ggrs = -2.d0*a*a1*q2-q0*q3/(q1**2+q1)
      return
      end
!
!  gga91 correlation
!  input rs: seitz radius
!  input zet: relative spin polarization
!  input t: abs(grad d)/(d*2.*ks*g)
!  input uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*g)**3)
!  input vv: (laplacian d)/(d * (2*ks*g)**2)
!  input ww: (grad d)*(grad zet)/(d * (2*ks*g)**2
!  output h: nonlocal part of correlation energy per electron
!  output dvcup,dvcdn:  nonlocal parts of correlation potentials
!

!*****************************************************************

      subroutine corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn,
     1                  fk,sk,g,ec,ecrs,eczet)
!
      implicit none
      real*8 rs,zet,t,uu,vv,ww,h,dvcup,dvcdn
      real*8 xnu,cc0,cx,alf,c1,c2,c3,c4,c5,c6,a4,thrdm,thrd2,bet,delt
      real*8 g3,g4,B,B2,BEC,BG,CC,CCRS,COEFF,COMM
      real*8 FAC,FACT0,FACT1,FACT2,FACT3,FACT4,FACT5
      real*8 GZ,H0,H0B,H0BT,H0RS,H0RST,H0T,H0TT,H0Z,H0ZT,H1,H1RS
      real*8 H1RST,H1T,H1TT,H1Z,H1ZT,HRS,HRST,HT,HTT,HZ,HZT
      real*8 PON,PREF,Q4,Q5,Q6,Q7,Q8,Q9,R0,R1,R2,R3,R4
      real*8 RS2,RS3,RSTHRD,T2,T4,T6
      real*8 fk,sk,g,ec,ecrs,eczet

      data xnu,cc0,cx,alf/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
      data c1,c2,c3,c4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
      data c5,c6,a4/0.472d0,7.389d-2,100.d0/
      data thrdm,thrd2/-0.333333333333d0,0.666666666667d0/
      bet = xnu*cc0
      delt = 2.d0*alf/bet
      g3 = g**3
      g4 = g3*g
      pon = -delt*ec/(g3*bet)
      b = delt/(dexp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      t6 = t4*t2
      rs2 = rs*rs
      rs3 = rs2*rs
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      q6 = c1+c2*rs+c3*rs2
      q7 = 1.d0+c4*rs+c5*rs2+c6*rs3
      cc = -cx + q6/q7
      r0 = (sk/fk)**2
      r1 = a4*r0*g4
      coeff = cc-cc0-3.d0*cx/7.d0
      r2 = xnu*coeff*g3
      r3 = dexp(-r1*t2)
      h0 = g3*(bet/delt)*log(1.d0+delt*q4*t2/q5)
      h1 = r3*r2*t2
      h = h0+h1
!  local correlation option:
!     h = 0.0d0
!  energy done. now the potential:
      ccrs = (c2+2*c3*rs)/q7 - q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
      rsthrd = rs/3.d0
      r4 = rsthrd*ccrs/coeff
!
      if(abs(zet) .ge. 1.d0) then
        if(zet .lt. 0.d0) zet=-1.d0+1.d-15
        if(zet .gt. 0.d0) zet=1.d0-1.d-15
      endif
      gz = ((1.d0+zet)**thrdm - (1.d0-zet)**thrdm)/3.d0
      fac = delt/b+1.d0
      bg = -3.d0*b2*ec*fac/(bet*g4)
      bec = b2*fac/(bet*g3)
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      h0b = -bet*g3*b*t6*(2.d0+b*t2)/q8
      h0rs = -rsthrd*h0b*bec*ecrs
      fact0 = 2.d0*delt-6.d0*b
      fact1 = q5*q9+q4*q9*q9
      h0bt = 2.d0*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
      h0rst = rsthrd*t2*h0bt*bec*ecrs
      h0z = 3.d0*gz*h0/g + h0b*(bg*gz+bec*eczet)
      h0t = 2*bet*g3*q9/q8
      h0zt = 3.d0*gz*h0t/g+h0bt*(bg*gz+bec*eczet)
      fact2 = q4*q5+b*t2*(q4*q9+q5)
      fact3 = 2.d0*b*q5*q9+delt*fact2
      h0tt = 4.d0*bet*g3*t*(2.d0*b/q8-(q9*fact3/q8)/q8)
      h1rs = r3*r2*t2*(-r4+r1*t2/3.d0)
      fact4 = 2.d0-r1*t2
      h1rst = r3*r2*t2*(2.d0*r4*(1.d0-r1*t2)-thrd2*r1*t2*fact4)
      h1z = gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
      h1t = 2.d0*r3*r2*(1.d0-r1*t2)
      h1zt = 2.d0*gz*r3*r2*(3.d0-11.d0*r1*t2+4.d0*r1*r1*t4)/g
      h1tt = 4.d0*r3*r2*r1*t*(-2.d0+r1*t2)
      hrs = h0rs+h1rs
      hrst = h0rst+h1rst
      ht = h0t+h1t
      htt = h0tt+h1tt
      hz = h0z+h1z
      hzt = h0zt+h1zt
      comm = h+hrs+hrst+t2*ht/6.d0+7.d0*t2*t*htt/6.d0
      pref = hz-gz*t2*ht/g
      fact5 = gz*(2.d0*ht+t*htt)/g
      comm = comm-pref*zet-uu*htt-vv*ht-ww*(hzt-fact5)
      dvcup = comm + pref
      dvcdn = comm - pref
!  local correlation option:
!     dvcup = 0.0d0
!     dvcdn = 0.0d0
      return
      end
