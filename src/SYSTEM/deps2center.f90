! deps2center.f90
! ===========================================================================
! This subroutine sets up deps/dr1 in the two-center molecular system 
! defined by sighat=(r2-r1)/|r2-r1| and piprimehat=(r2-cross-r1)/|r2-cross-r1|.
! The eps matrix is obtained by calling epsiln(r1,rvec) where rvec=r2-r1.
! ===========================================================================
subroutine deps2cent(r1,r2,eps2,deps2)
  use M_constants_fireball
  implicit none

  ! input:
  ! r1    - position of atom 1
  ! r2    - posiiton of atom 2
  ! eps2  - 2 center epsilon with z=r2-r1
  real, intent(in) :: r1(3),r2(3),eps2(3,3)

  ! output:
  ! deps2 - deps/dr1
  real, intent(out) :: deps2(3,3,3)
 
  integer i,ii,ix

  real r2mag2,r2mag,r1mag,denom
  real crossmag,dd,dot,term,ddinv,crossinv
  real crossa(3)
 
  deps2=0.e0

  ! If we are doing an atom, r1=r2. Then set deps2 to zero.
  dd=sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2) 
  if(dd.lt.1.0d-4)return
  ddinv=1.0/dd 
 
  r2mag2= r2(1)**2 + r2(2)**2 + r2(3)**2
  r2mag=sqrt(r2mag2)
  r1mag=sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
  crossa(1)=r2(2)*r1(3)-r2(3)*r1(2)
  crossa(2)=r2(3)*r1(1)-r2(1)*r1(3)
  crossa(3)=r2(1)*r1(2)-r2(2)*r1(1)
  crossmag=sqrt(crossa(1)**2+crossa(2)**2+crossa(3)**2)
  dot=0.e0
  do i=1,3
    dot=dot+r1(i)*r2(i)
  end do
  denom=r1mag*r2mag
 
  ! check to see if any atom is near origin
  if(denom.lt.1.0d-3)then
    return
  endif
  ! check to see if atoms are colinear
  if(abs(crossmag).lt.1.0d-3)then
    return
  endif
  crossinv=1.0/crossmag
  ! now calculate deps
  do ii=1,3
   do ix=1,3
    term=xlevi(ix,ii,1)*r2(1)+ xlevi(ix,ii,2)*r2(2)+xlevi(ix,ii,3)*r2(3)

    deps2(ix,ii,1)=(eps2(ii,1)*eps2(ix,3)*ddinv)-(eps2(ii,1)*crossinv**2)*(r2mag2*r1(ix)-dot*r2(ix))+      &
     &              ddinv*crossinv*(delk(ii,ix)*(r2mag2-dot)-r1(ii)*r2(ix)-r2(ii)*r2(ix)+2.e0*r1(ix)*r2(ii))

    deps2(ix,ii,2)=crossinv*(term+(eps2(ii,2)*crossinv)*(dot*r2(ix)-r2mag2*r1(ix)))
    deps2(ix,ii,3)=-(delk(ii,ix)-eps2(ii,3)*eps2(ix,3))*ddinv
   end do
  end do
  return
end
