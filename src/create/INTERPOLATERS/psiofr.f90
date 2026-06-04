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
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.


! psiofr.f
! Program Description
! ===========================================================================
!       This function returns the values psiofr(r) for the corresponding
! shell of the atomtype ispec.  The radial functions are normalized as:

!  int ( psiofr**2  r**2  dr ) = 1.0

! The wavefunctions for each atom must be read in by calling readpsi
! separately for each atom type.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function psiofr (ispec, issh, r)
        use precision, only: wp
        implicit none
        real(kind=wp) psiofr

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
        integer ispec
        integer issh

        real(kind=wp) r

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

        integer norder2
        parameter (norder2 = norder/2)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer max_points

        real(kind=wp) L(0:norder),mu(0:norder),Z(0:norder),alpha(0:norder)
        real(kind=wp) a(0:norder),b(0:norder),c(0:norder),d(0:norder)
        real(kind=wp) h
        integer iam
        integer i,j
        real(kind=wp) xmin
        real(kind=wp) xxp
        real(kind=wp) aaa,bbb,ccc,ddd

! Procedure
! ===========================================================================
! Special cases
        if (r .ge. rrc(issh,ispec)) then
          psiofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          psiofr = psi(1,issh,ispec)
          return
        end if

! note : the points are equally spaced
        h=drr(issh,ispec)

        imid = int(r/h) + 1
        max_points = npoints(issh,ispec)

        if (superspline) goto 1234

! Find starting point for the interpolation
        ileft = imid - norder2
        if (ileft .lt. 1) then
         ileft = 1
        else if (ileft + norder .gt. max_points) then
         ileft = max_points - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points

        do i=0,norder
          a(i)=psi(i+ileft,issh,ispec)
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

!       What curve section do we use?
        iam=imid-ileft

!       Don't need 0 to iam-1
        do j=norder-1,iam,-1
          c(j)=z(j)-mu(j)*c(j+1)
          b(j)=(a(j+1)-a(j))/h-h*(c(j+1)+2.0d0*c(j))/3.0d0
          d(j)=(c(j+1)-c(j))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp=(r-(xmin+(imid-1)*h))

        psiofr=a(iam)+b(iam)*xxp+c(iam)*xxp**2+d(iam)*xxp**3

        return
 
!
! Cubic splines:  One big "super spline"
!
 1234   continue 
        if (imid .gt. max_points) imid=max_points ! If we have gone off of the end
        xxp=r-(imid-1)*h
        psiofr = psi_spline(4,imid,issh,ispec)
        do j = 3, 1, -1
          psiofr = xxp*psiofr + psi_spline(j,imid,issh,ispec)
        end do

! Format Statements
! ===========================================================================

        return
        end
