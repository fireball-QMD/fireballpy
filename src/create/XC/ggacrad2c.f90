!**********************************************************************
!
! correlation potential & energy density
! spherical symmetry
! LSDA - GGA
! Hartree a.u.
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-C Perdew 91
!    mode = 3    GGA-C Perdew 86
!    mode = 4    GGA-C Lee-Yang-Parr 1988
!    mode = 5    GGA-C Burke-Perdew-Ernzerhof
!    r           radius
!    rh()        spin up/down density
!    rhp()       1st derivative of rh
!    rhpp()      2nd derivative of rh
!
! output
!    zet         spin polarization
!    cpot()       - " -  correlation potential
!    cen         correlation energy density
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!**********************************************************************
!
      subroutine ggacrad2c (mode, r, rh, rhp, rhpp, rhz, rhzz, &
     &                      cpot, cen)
!
      use precision
      implicit none

      real(kind=long) r, rh, rhp, rhpp, rhz, rhzz, cpot, cen
      real(kind=long) thrd,pi,pisq3,thrd2,crs,eps
      real(kind=long) alfc,d,dp,dp11,dp12,dp22,dpp,dvcdn,dvcup,gks2
      real(kind=long) gks2sq,h,rs,t,uu,vcdn,vcup,vv,ww,zet,ztp
      real(kind=long) fk,sk,g,ec,ecrs,eczet

      integer mode

      dimension rh(2),rhp(2),rhpp(2),cpot(2)
      dimension rhz(2),rhzz(2)
      data thrd,pi /.333333333333333333d0,.314159265358979312d1/
      data pisq3,thrd2 /.296088132032680740d2,.666666666666666666d0/
      data crs,eps /1.91915829267751281d0,1.d-15/
!
! LSDA
      d = rh(1)+rh(2)
      cen=0.d0
      if(d .le. eps) then
        cen=0.d0
        cpot(1)=0.d0
        cpot(2)=0.d0
      else
        if(mode .ne. 4) then
          zet=(rh(1)-rh(2))/d
          fk=(pisq3*d)**thrd
          rs=crs/fk
!
! GGA correction to LSDA
          if(mode .eq. 1) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dvcup=0.d0
            dvcdn=0.d0
            h=0.d0
          else if(mode .eq. 2) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dp=rhp(1)+rhp(2)
            dpp=rhpp(1)+rhpp(2)
            ztp=(rhp(1)-rhp(2)-zet*dp)/d
            sk=2.d0*sqrt(fk/pi)
            g=((1.d0+zet)**thrd2 + (1.d0-zet)**thrd2) /2.d0
            gks2=2.d0*sk*g
            gks2sq=gks2*gks2
            t=abs(dp)/(d*gks2)
            uu=abs(dp)*dpp/(d*d*gks2sq*gks2)
            vv=(dpp +2.d0*dp/r)/(d*gks2sq)
            ww=dp*ztp/(d*gks2sq)
            call corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn, &
     &                  fk,sk,g,ec,ecrs,eczet)
          else if(mode .eq. 3) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dp  =rhp(1)+rhp(2)
            dpp =rhpp(1)+rhpp(2)
            uu  =abs(dp)*dpp
            vv  =dpp +2.d0*dp/r
            dp11=rhp(1)*rhp(1)
            dp22=rhp(2)*rhp(2)
            dp12=rhp(1)*rhp(2)
            if(rhp(1) .ne. 0.d0 .or. rhp(2) .ne. 0.d0) &
     &      call corga86(rh(1),rh(2),dp11,dp22,dp12,uu,vv,h,dvcup,dvcdn)
          else if(mode .eq. 5) then
            dp=rhp(1)+rhp(2)
            dpp=rhpp(1)+rhpp(2)
            ztp=(rhp(1)-rhp(2)-zet*dp)/d
            sk=2.d0*sqrt(fk/pi)
            g=((1.d0+zet)**thrd2 + (1.d0-zet)**thrd2) /2.d0
            gks2=2.d0*sk*g
            gks2sq=gks2*gks2
            t=abs(dp)/(d*gks2)
            uu=abs(dp)*dpp/(d*d*gks2sq*gks2)
            vv=(dpp +2.d0*dp/r)/(d*gks2sq)
            ww=dp*ztp/(d*gks2sq)
            call CORPBE(rs,zet,t,uu,vv,ww,1,1, &
     &                     ec,vcup,vcdn,h,dvcup,dvcdn)
          else
            stop 'ggacrad : mode improper'
          endif
          cpot(1)=vcup+dvcup
          cpot(2)=vcdn+dvcdn
          cen=ec+h
        else if(mode .eq. 4) then
          call corlyp2c (.true., r, rh(1), rh(2), rhp(1), rhp(2), &
     &                   rhpp(1), rhpp(2), rhz(1), rhz(2), rhzz(1), &
     &                   rhzz(2), cen, cpot(1))
        else
          stop 'ggacrad : mode improper'
        endif
      endif
!
      return
      end
!
