! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!                      Christopher Sikorski
!                      Xiaodong Chen
!
! Universidad Autonoma de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
!
! Code changed by Jose Ortega : the 3-center integral is now performed
!                               with the origin at the center of the
!                                          bond-charge 

! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
!     ============================================================
       subroutine threecenter_integral(dbcx,rna,rc1,rc2,nrr,nt,nphi,
     1                   in1,in2,in3,gmat,index_max,n1,l1,m1,n2,l2,m2,
     2                   iexc,interaction,ispmin,ispmax,ispher)
!     ============================================================
!
!     Calculates the matrix elements of the neutral atom for the
!     specific configuration of dbcx (dist. between bc) and rna=position of
!     neutral atom.
!
!     ------------------------------------------------------------
!
      implicit none
!
      include '../parameters.inc'
      include '../vnonneutral.inc'
      include '../wavefunctions.inc'

!
!     ------------------------------------------------------------
!
!     Passed variables:
      integer numInthpR, numinthpT, numinthpP
!
! JOM : the following parameters are the grids for integration
! JOM : I think these numbers should be of the
!       type (6*N +1)  (Newton-Cotes 7-point rule) 
! so far I recommend (note that symmetry (y,-y) is now used
! in the phi-integral
      parameter (numInthpR=49)
      parameter (numInthpT=49)
      parameter (numInthpP=25)
!
      integer in1,in2,in3      ! types of the atoms
      integer iexc             ! exchange-correlation approximation
      integer nrr,nt,nphi      ! number of integration points
!
      integer index_max        ! maximal number of matrix elements
      integer interaction      ! number for the interaction
      integer ispmin,ispmax    ! range for isorp loop
      logical ispher           ! spherical approx. 
!
      integer n1 (inter_max)   ! shell number of left  atom
      integer n2 (inter_max)   !                 right atom
      integer l1 (inter_max)   ! angular momentum of that shell
      integer l2 (inter_max)   !
      integer m1 (inter_max)   ! m-value in that shell
      integer m2 (inter_max)   !
!
      real*8 rc1           ! i=1,2,3  for left,right, neutral atm
      real*8 rc2
!
      real*8 h
      real*8 HR,HT,HP
      integer tempR, tempT, tempP
 
      real*8 dbcx          ! bond charge distance  (A)
      real*8 rna(3)        ! neutral atom location (A)
!
      integer tmpImid,imid
      real*8 xmin
      real*8 xxp
      real*8 psipsi
!
      real*8 gmat(0:10,inter_max)   ! result from the integrator
!
      real*8 psiofr,vpot,vnnaofr,dvxc3c
      external dvxc3c, vnnaofr,psiofr
!
!     ............................................................
!
!     interaction = 1: bcna  neutral atom
!                   2: xc3c  exchange correlation
!                   3: xc3c  exchange correlation (SNXC and OLSXC)
!
!     for the exchange correaltion case,
!
!     gmat(ix,*) stores the derivatives w.r.t charges:
!
!     gmat(1,*) : in1,0    gmat(2,*) : in1,-q1  gmat(3,*) : in1,+q1
!     gmat(4,*) : in2,-q2  gmat(5,*) : in2,+q2
!     gmat(6,*) : in3,-q3  gmat(7,*) : in3,+q3
!
!     ............................................................
!
!     ------------------------------------------------------------
!
!     Internal variables:
!
!
      real*8 thetamat(0:10,inter_max)
!
      real*8 znormMU(inter_max),znormNU(inter_max),
     1       psiAmat(inter_max),psiBmat(inter_max)
!
!
      real*8 thfactMU(0:3,-3:3),
     1       thfactNU(0:3,-3:3),
     2       phifactor(-3:3)
!
      real*8 avgVmat(0:10,inter_max), fsimp(5000)
!
      real*8 inthpR(numInthpR*2+1), inthpTheta(numInthpT*2+1), 
     1 inthpPhi(numInthpP*2+1)
      real*8 wR(numInthpR*2+1),wTheta(numInthpT*2+1),
     1 wPhi(numInthpP*2+1)
!
      integer   numbphi
      integer irMax, itMax,iPhiMax
      parameter (numbphi=5000)
      real*8  paramR, paramT, paramP
      real*8    phiy(numbphi),
     1          cphiy(numbphi),sphiy(numbphi)

      real*8 w1, w2
!
      integer nn1,nl1,nm1,nn2,nl2,nm2,i,ir,ix,ix1,it,ip,inm,
     1        nmax,isorpX, NinthpR, NinthpT, NinthpP
!
! JOM r1
      real*8  r1
      real*8  sq3,sq15,pi,dr,dtheta,dphi,r,theta,phi,
     1        rmin,rmax,dc1,ds1,dc,ds,zr,r2,dc2,ds2,xr,yr,r3,
     2        cphi,sphi,simp,simp2,simpson,
     3        averagephi,averagetheta,
     4        prod,prod2,dsth,stuffmunu,thrd,nrrinv,ntinv,nphiinv
      real*8 sq12,sq4,hunderedfortieth
!
! =====================================================================
!     The size of the matrix is determined by Nsh(in1) and Nsh(in2)
!     We do only those matrix elements that are not zero
!     The number of the non-zero matrix elements is INMAX
!
!     Normalization factors
      
      sq3=1.73205080756887729352744634150587d0
      sq15=3.87298334620741688517926539978240d0
      sq12=3.46410161513775458705489268301174d0
      sq4 = 2.0d0
      pi= 3.14159265358979323846264338327950d0
      thrd=0.333333333333333333333333333333333d0

      rmin=0.0d0

      ! added my murat manguoglu
      hunderedfortieth=0.007142857142857142857142857142857143d0 ! 1/140

! JOM
!      rmax = dmax1(rc1,rc2)
      rmax = 0.5d0*(rc1+rc2)

      NinthpR = numInthpR
      NinthpT = numInthpT
      NinthpP = numInthpP
    
      !modified by murat manguoglu 
      irMax = NinthpR
      itMax = NinthpT
      iphiMax = NinthpP
      ! end of modification by murat manguoglu 
      
      paramR = 2.0D0  ! 2*d/pi
      paramT = 2.0D0
      paramP = 1.0D0
  
      nrrinv=1.0d0/dfloat(irMax-1)
      ntinv=1.0d0/dfloat(itMax-1)
      nphiinv=1.0d0/dfloat(iphiMax-1)

      HR = (rmax-rmin)*nrrinv
      HT = PI*ntinv
! JOM-test
!      HP = 2.D0*PI*nphiinv
      HP = PI*nphiinv
! JOM
! I think I will change the integral in phi from (0,2pi) to (0,pi): latter

! now we choose d = Pi/4      
      
! modified by murat manguoglu
  
      inthpR(1)=rmin
      wR(1) = HR*hunderedfortieth*41.0D0 
      do ir= 2, irMax-1
         inthpR(ir) = rmin+(ir-1.0D0)*HR
           if (mod(ir,6).eq.2) wR(ir)=HR*hunderedfortieth*216.0D0 
           if (mod(ir,6).eq.3) wR(ir)=HR*hunderedfortieth*27.0D0 
           if (mod(ir,6).eq.4) wR(ir)=HR*hunderedfortieth*272.0D0
           if (mod(ir,6).eq.5) wR(ir)=HR*hunderedfortieth*27.0D0
           if (mod(ir,6).eq.0) wR(ir)=HR*hunderedfortieth*216.0D0
           if (mod(ir,6).eq.1) wR(ir)=HR*hunderedfortieth*82.0D0 
      end do
      wR(irMax)= HR*hunderedfortieth*41.0D0
      inthpR(irMax)= rmax
    

      inthpTheta(1)=0.0D0
      wTheta(1) = HT*hunderedfortieth*41.0D0
      do it= 2, itMax-1
         inthpTheta(it) = (it-1.0D0)*HT
           if (mod(it,6).eq.2) wTheta(it)=HT*hunderedfortieth*216.0D0
           if (mod(it,6).eq.3) wTheta(it)=HT*hunderedfortieth*27.0D0
           if (mod(it,6).eq.4) wTheta(it)=HT*hunderedfortieth*272.0D0
           if (mod(it,6).eq.5) wTheta(it)=HT*hunderedfortieth*27.0D0
           if (mod(it,6).eq.0) wTheta(it)=HT*hunderedfortieth*216.0D0
           if (mod(it,6).eq.1) wTheta(it)=HT*hunderedfortieth*82.0D0
      end do
      wTheta(itMax)= HT*hunderedfortieth*41.0D0
      inthpTheta(itMax)= PI


      inthpPhi(1)=0.0D0
      wPhi(1) = HP*hunderedfortieth*41.0D0
      do ip= 2, iphiMax-1
         inthpPhi(ip) = (ip-1.0D0)*HP
         if (mod(ip,6).eq.2) wPhi(ip)=HP*hunderedfortieth*216.0D0
         if (mod(ip,6).eq.3) wPhi(ip)=HP*hunderedfortieth*27.0D0
         if (mod(ip,6).eq.4) wPhi(ip)=HP*hunderedfortieth*272.0D0
         if (mod(ip,6).eq.5) wPhi(ip)=HP*hunderedfortieth*27.0D0
         if (mod(ip,6).eq.0) wPhi(ip)=HP*hunderedfortieth*216.0D0
         if (mod(ip,6).eq.1) wPhi(ip)=HP*hunderedfortieth*82.0D0
      end do
      wPhi(iphiMax)= HP*hunderedfortieth*41.0D0
! JOM-test
!      inthpPhi(iphiMax)= 2.0D0*PI
      inthpPhi(iphiMax)= PI

!end of modification by murat manguoglu

! ===================================================================
! Here is the correct list:
!   m:     -2       -1         0        1         2
!          xy       yz      3z^2-r^2   xz       x^2-y^2
! sq15 * (  1        1        1/sq12    1       1/sq4  )
!
! ====================================================================
      do 313 inm=1,index_max
        nl1=l1(inm)
        nl2=l2(inm)
        if(nl1.eq.0)znormMU(inm)=1.0d0
        if(nl2.eq.0)znormNU(inm)=1.0d0
        if(nl1.eq.1)znormMU(inm)=sq3
        if(nl2.eq.1)znormNU(inm)=sq3
        if(nl1.eq.2)znormMU(inm)=sq15
        if(nl2.eq.2)znormNU(inm)=sq15
! First we calculate the m values.
        nm1=m1(inm)
        nm2=m2(inm)
! working on d for mu
        if(nl1.eq.2)then
          if(nm1.eq.0)znormMU(inm)=znormMU(inm)/sq12
          if(nm1.eq.2)znormMU(inm)=znormMU(inm)/sq4
        end if
! working on d for nu
        if(nl2.eq.2)then
          if(nm2.eq.0)znormNU(inm)=znormNU(inm)/sq12
          if(nm2.eq.2)znormNU(inm)=znormNU(inm)/sq4
        end if
 313  continue
 
!
!     isorp=0,1,2,3, for neutral atom, s part, p part, d part.
!     We add 1 because of 0 being the neutral atom.
!
 
      if(nphi.gt.numbphi)then
        write(*,*)' nphi=',nphi
        write(*,*)' numbphi=',numbphi
        write(*,*)' In threecenter_integral-- redimension numbphi'
        stop 'error in threecenter_integral'
      end if
!
! ========================================================
!

! JOM we do not use the following lines:
!      nmax=max(nrr,nt,nphi)
!
! ========================================================
!
!     rmax=rc1
!      rmax = dmax1(rc1,rc2)
!
!      dr=(rmax-rmin)*nrrinv
!      dtheta=pi*ntinv
!      dphi=2.d0*pi*nphiinv
!
! end JOM

!
!     Set up some constants.
!     The phi integration does not depend on the r integration.
!     Therefore, it is done outside the r loop.
!
      do 7244 ip= 1, iphiMax
        phi=inthpPhi(ip)
        cphiy(ip)=cos(phi)
        sphiy(ip)=sin(phi)
 7244 continue
! ========================================================
 
      do 1451 inm=1,index_max
        do 1451 isorpX=ispmin,ispmax
          gmat(isorpX,inm)=0.0d0
 1451 continue
     
!
! ========================================================
!                 Do integral over r:
! ========================================================
!

 
        do 310 ir= 1, irMax
        r=inthpR(ir)
!
!        dorbital addition OFS
!        Fireball2000: Jose y Jandro.
!        Zero out the array for theta integral.
!
         do inm=1,index_max
           do ix=ispmin,ispmax
             thetamat(ix,inm)=0.0d0
           end do
         end do
!
!        nli, nni tells us which shell is used for each atom.
!

! JOM we have to define psiAmat latter
!         do inm=1,index_max
!           nn1=n1(inm)
!           psiAmat(inm)=psiofr(in1,nn1,r)
!! jel-spher
!           if(ispher) psiAmat(inm) = sqrt(psiAmat(inm)**2.0d0)
!         end do
! end JOM

 
!
! ========================================================
!              Do integral over theta:
! ========================================================
!
         do 208 it=1, itMax 
!           
           theta=inthpTheta(it)
!          Theta stuff for atom 1

!=========================================================
! JOM JOM  : main changes
!=========================================================
!           dc1=cos(theta)
!           ds1=sin(theta)
!           dc=dc1
!           ds=ds1
           dc=cos(theta)
           ds=sin(theta)
!
	   r1 = r**2 + 0.25d0*(dbcx**2) + r*dbcx*dc
           if (r1 .le. 0.0d0) then
            r1 = 0.0d0
           else
            r1 = sqrt(r1)
           end if
           if(r1.gt.rc1)go to 208       ! outside integration range

!           r2 = r**2 + dbcx**2 - 2*r*dbcx*dc
	   r2 = r**2 + 0.25d0*(dbcx**2) - r*dbcx*dc
           if (r2 .le. 0.0d0) then
            r2 = 0.0d0
           else
            r2 = sqrt(r2)
           end if
           if(r2.gt.rc2)go to 208       ! outside integration range

!                                       ! other theta values might work
!
           do inm=1,index_max
             nn1 = n1(inm)
             nn2 = n2(inm)
             psiAmat(inm)=psiofr(in1,nn1,r1)
             psiBmat(inm)=psiofr(in2,nn2,r2)
! jel-spher
             if(ispher) then
              psiAmat(inm) = sqrt(psiAmat(inm)**2.0d0)
              psiBmat(inm) = sqrt(psiBmat(inm)**2.0d0)
              end if
           enddo
!
           zr=r*dc
!          Theta stuff for atom 1.
!          Be careful for r1 very small.
!          Find cos(theta1), sin(theta1).
           if(r1.gt.0.00001)then
             dc1=(zr+0.5d0*dbcx)/r1
             ds1=ds*r/r1
           else
             dc1=1.0d0
             ds1=0.0d0
           end if
!
!          Theta stuff for atom 2.
!          Be careful for r2 very small.
!          Find cos(theta2), sin(theta2).
           if(r2.gt.0.00001)then
             dc2=(zr-0.5d0*dbcx)/r2
             ds2=ds*r/r2
           else
             dc2=1.0d0
             ds2=0.0d0
           end if
!=========================================================
! end JOM JOM
!=========================================================
!
!
! -------------------------------------------------------
!          Theta factors for A and B.
! -------------------------------------------------------
!          Use the (l,m) notation. For example 3z^2-1 becomes
!          (2,0), px becomes (1,1), and so on.
!
!          ATOM A .........................................
!
!          S
           thfactMU(0,0)=1.0d0
!
!          P
!          Note: We order the orbitals here x,y,z (or pi,pi',sig)
!
!          watch it: theta factors for d-orbitals should
!                    also contain the m-dependency of the
!                    prefactors znormMU and znormNU.
!                    This is not the case so far.
!
           thfactMU(1,1)=ds1
           thfactMU(1,-1)=ds1
           thfactMU(1,0)=dc1
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactMU(2,0)=3.0d0*dc1*dc1-1.0d0
           thfactMU(2,2)=ds1*ds1
           thfactMU(2,1)=ds1*dc1
           thfactMU(2,-2)=ds1*ds1
           thfactMU(2,-1)=ds1*dc1
!
!          ATOM B .............................................
!
!          S
           thfactNU(0,0)=1.0d0
!
!          P
!
           thfactNU(1,1)=ds2
           thfactNU(1,-1)=ds2
           thfactNU(1,0)=dc2
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactNU(2,0)=3.0d0*dc2*dc2-1.0d0
           thfactNU(2,2)=ds2*ds2
           thfactNU(2,1)=ds2*dc2
           thfactNU(2,-2)=ds2*ds2
           thfactNU(2,-1)=ds2*dc2
!
! --------------------------------------------------------------
!          Done with theta factors.
! --------------------------------------------------------------
 
           do inm=1,index_max
             do ix=ispmin,ispmax
               avgVmat(ix,inm)=0.0d0
             end do
           end do
!
! ========================================================
!          Do integral over phi:
! ========================================================
 
!          average over phi (divide by 2*pi)
! JOM : I add this normalization factor at the end, out
!       of the loops
!           averagephi=0.5d0/pi
! HAO : It should be averagephi = 1.d0/pi, since interal of phi is now
!       done in (0, pi).
 
!          The phi factors depend only on m.
 
           phifactor(0)=1.0d0
 
           do 244 ip=1, iphiMax
!

             phi=inthpPhi(ip)
             cphi=cphiy(ip)
             sphi=sphiy(ip)
!
!            Do the phi integral, with the phifactors.
!            Note: We order the p-orbitals
!            here x,y,z (or pi,pi',sig), NOT z,x,y.
!            Note that px, and xz now are +1. And so on!
! JOM : I think the above is not true: y,z,x (-1,0,1)
!
             phifactor(1)=cphi
             phifactor(-1)=sphi
! d's
             phifactor(2)=cphi*cphi-sphi*sphi
             phifactor(-2)=cphi*sphi
 
             xr=r*ds*cphi
             yr=r*ds*sphi
!
             r3=sqrt((xr-rna(1))**2+(yr-rna(2))**2+(zr-rna(3))**2)
!
! ---------------------------------------------------
!
             do  iX=ispmin,ispmax
!
               IF(interaction .EQ. 1) vpot=vnnaofr(in3,iX,r3)
!
               IF(interaction .EQ. 2) then
                IX1=IX+1
                vpot = dvxc3c (iexc, r1, r2, r3, in1, in2, in3, IX1)
               END IF
! xc3c_SN
               IF(interaction .EQ. 3) then
                  psipsi = psiofr(in3,ix,r3)
                  vpot=(psipsi**2)/(4.0d0*pi)
               ENDIF
!
!              note: dc,ds defined at the beginning of the theta loop
!

! JOM : I add this normalization factor at the end, out
!       of the loops
!               prod=vpot*wPhi(ip)*averagephi
               prod=vpot*wPhi(ip)
               do inm=1,index_max
                 nm1=m1(inm)
                 nm2=m2(inm)
                 avgVmat(ix,inm)=avgVmat(ix,inm) +
     1                           prod*phifactor(nm1)*phifactor(nm2)
               end do
             end do
!
! ===============================================
!          The end of the phi integral.
 244       continue
! ===============================================
!
! JOM
!           dsth=ds1
           dsth=ds
! JOM I add this norm. factor outside loops
!           averagetheta=0.5d0
!           prod=wTheta(it)*averagetheta*dsth
           prod=wTheta(it)*dsth
 
           do 1407 inm=1,index_max
             nl1=l1(inm)
             nm1=m1(inm)
             nl2=l2(inm)
             nm2=m2(inm)
! JOM add psiAmat here
!             stuffmunu=prod*thfactMU(nl1,nm1)*
!     &                      thfactNU(nl2,nm2)*psiBmat(inm)
             stuffmunu=prod*thfactMU(nl1,nm1)*
     1              thfactNU(nl2,nm2)*psiBmat(inm)*psiAmat(inm)
             do 1408 ix=ispmin,ispmax
               thetamat(ix,inm)=thetamat(ix,inm)+
     1                          avgVmat(ix,inm)*stuffmunu
 1408        continue
 1407      continue
!
!
!
 208    continue
!       ========================================
!       The end of the integral over theta (loop 208).
!       =====================================================
 
!       now finish off the r integral!
 
        prod=wR(ir)*r*r
 
        do inm=1,index_max
! JOM remove psiAmat form here
!          prod2=prod*psiAmat(inm)
          do ix=ispmin,ispmax
! JOM
!            gmat(ix,inm)=gmat(ix,inm)+prod2*thetamat(ix,inm)
            gmat(ix,inm)=gmat(ix,inm)+prod*thetamat(ix,inm)
          end do
        end do
!
!
!
 310  continue
!
! ========================================================
!     The end of the integral over r (loop 310).
! ========================================================
!
!     Finally, the normalization factors for the different orbitals.
!     For instance a p orbital is sqrt(3) * x/r * r(r). The sqrt(3)
!     factor (and ones like it) are now included.
! 
! JOM also averagetheta*averagephi = 1/(2 pi) is added now
!
      do inm=1,index_max
!        prod=znormMU(inm)*znormNU(inm)
!JOM
        prod=znormMU(inm)*znormNU(inm)*0.5d0/pi
        do ix=ispmin,ispmax
           gmat(ix,inm)=gmat(ix,inm)*prod
        end do
      end do
!
! ================================================================
!     SUMMARY
! We have computed gmat(isorp,mu,nu). mu, and nu are 1 to 9 for
! sp^3d^5.
! We are in molecular coordinates with z along sigma, x along pi, and
! y along pi'. Many matrix elements og gmat are zero. We have computed them
! anyway, and indeed find they are zero. Just to avoid at a later time,
! any trouble with roundoffs, I will now set those that are supposed to be
! zero, identically to zero.
! Here is the matrix, and the zero's are indicated.
!
!          \ s   x   y   z    3z^2-r^2    x^2-y^2    xz    xy      yz
!
!  s                 0                                     0        0
!
!  x                 0                                     0        0
!
!  y         0   0       0        0          0        0
!
!  z                 0                                     0        0
!
! 3z^2-r^2           0                                     0        0
!
! x^2-y^2            0                                     0        0
!
!  xz                0                                     0        0
!
!  xy        0   0       0        0          0        0
!
!  yz        0   0       0        0          0        0
!
!
       return
       end


