subroutine deps2cent(r1,r2,eps2,deps2)
  use M_constants, only: wp, xlevi, delk

  implicit none
  real(wp), dimension(3), intent(in) :: r1, r2
  real(wp), dimension(3, 3), intent(in) :: eps2
  real(wp), dimension (3, 3, 3), intent(out) :: deps2
  integer :: i,ii,ix
  real(wp) :: r2mag2,r2mag,r1mag,denom,crossmag,dd,dot,term,ddinv,crossinv
  real(wp), dimension(3) :: crossa
  deps2=0.e0
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
  if(denom.lt.1.0d-3)then
    return
  endif
  if(abs(crossmag).lt.1.0d-3)then
    return
  endif
  crossinv=1.0/crossmag
  do ii=1,3
   do ix=1,3
    term=xlevi(ix,ii,1)*r2(1)+ xlevi(ix,ii,2)*r2(2)+xlevi(ix,ii,3)*r2(3)
    deps2(ix,ii,1)=(eps2(ii,1)*eps2(ix,3)*ddinv)-(eps2(ii,1)*crossinv**2)*(r2mag2*r1(ix)-dot*r2(ix))+ ddinv*crossinv*(delk(ii,ix)*(r2mag2-dot)-r1(ii)*r2(ix)-r2(ii)*r2(ix)+2.e0*r1(ix)*r2(ii))
    deps2(ix,ii,2)=crossinv*(term+(eps2(ii,2)*crossinv)*(dot*r2(ix)-r2mag2*r1(ix)))
    deps2(ix,ii,3)=-(delk(ii,ix)-eps2(ii,3)*eps2(ix,3))*ddinv
   end do
  end do
  return
end
