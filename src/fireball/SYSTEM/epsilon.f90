! epsilon.f90
! Program Description
! ===========================================================================
! input: r1,r2
! output: epsilon the metric tensor
!
! note: the third column of epsilon is eta(3)
!
! spe=epsilon backwards.
!
! R1vector points toward the point O while R2vector points
!   away from the point O.
!
!                         *O    (XP,YP,ZP)
!                      *    *
!       R1VECTOR    *        *  R2VECTOR
!                *            *
!             *                *
! (X,Y,Z)  *                    *
!       *
!
!            |  ^     ^       ^     ^        ^     ^   |
!            |  X-dot-XP      X-dot-YP       X-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!      spe = |  Y-dot-XP      Y-dot-YP       Y-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!            |  Z-dot-XP      Z-dot-YP       Z-dot-ZP  |
!            |                                         |
!
! ===========================================================================
subroutine epsilon(R1,R2,spe)
  use iso_c_binding
  implicit none
 
  real(c_double), intent(in) :: r1(3)
  real(c_double), intent(in) :: r2(3)
  real(c_double), intent(out) :: spe(3,3)
 
  integer(c_long) i,j,ii,jj,kk,ix
  real(c_double) r1mag,r2mag,ypmag,unit
  real(c_double) XPHAT(3),YPHAT(3),ZPHAT(3),R1HAT(3)
 
  r1mag=sqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
  r2mag=sqrt(r2(1)*r2(1)+r2(2)*r2(2)+r2(3)*r2(3))
  if (r2mag .lt. 1.0d-4) then
    !r2vector = (0,0,0) ----- set eps = unit matrix
    do i=1,3
      do j=1,3
        spe(i,j)=0.0e0
      end do
      spe(i,i)=1.0e0
    end do
    return
  end if
  ! zphat lies along r2vector
  zphat(1)=r2(1)/r2mag
  zphat(2)=r2(2)/r2mag
  zphat(3)=r2(3)/r2mag

  ! yphat = zphat-cross-r1hat
  if(r1mag.gt.1.0d-4)then
    r1hat(1)=r1(1)/r1mag
    r1hat(2)=r1(2)/r1mag
    r1hat(3)=r1(3)/r1mag
    yphat(1)=zphat(2)*r1hat(3)-zphat(3)*r1hat(2)
    yphat(2)=zphat(3)*r1hat(1)-zphat(1)*r1hat(3)
    yphat(3)=zphat(1)*r1hat(2)-zphat(2)*r1hat(1)
    ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
    if(ypmag.gt.0.000001)goto 3000
  end if
  ! zphat and r1hat are colinear or r1vector=(0,0,0)
  !   find the first non-zero component of zphat
  if(abs(zphat(1)).gt.1.0d-4)then
  ! zphat(1) not equal to zero
    yphat(1)=-(zphat(2)+zphat(3))/zphat(1)
    yphat(2)=1.0e0
    yphat(3)=1.0e0
    ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
  else if(abs(zphat(2)).gt.1.0d-4)then
  ! zphat(2) not equal to zero
    yphat(1)=1.0e0
    yphat(2)=-(zphat(1)+zphat(3))/zphat(2)
    yphat(3)=1.0e0
    ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
  else
  ! zphat(3) not equal to zero
    yphat(1)=1.0e0
    yphat(2)=1.0e0
    yphat(3)=-(zphat(1)+zphat(2))/zphat(3)
    ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
  end if

  3000   continue

  yphat(1)=yphat(1)/ypmag
  yphat(2)=yphat(2)/ypmag
  yphat(3)=yphat(3)/ypmag
  ! find pihat
  xphat(1)=yphat(2)*zphat(3)-yphat(3)*zphat(2)
  xphat(2)=yphat(3)*zphat(1)-yphat(1)*zphat(3)
  xphat(3)=yphat(1)*zphat(2)-yphat(2)*zphat(1)
  ! find epsilon matrix
  do ix=1,3
    spe(ix,1)=xphat(ix)
    spe(ix,2)=yphat(ix)
    spe(ix,3)=zphat(ix)
  end do
 
  ! test by computing spe*spe(dagger)
  do ii=1,3
    do jj=1,3
      unit=0.0d0
      do kk=1,3
        unit=unit+spe(ii,kk)*spe(jj,kk)
      end do
    end do
  end do
  return
end subroutine epsilon
