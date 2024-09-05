subroutine deps3center (r1, r2, r21, distance12, ratm, rnabc, eps3, dera3, der13)
  use iso_c_binding
  use M_constants, only: xlevi, delk
  implicit none
  real(c_double), intent(in) :: distance12 ! distance between atom 1 and atom 2
  real(c_double), dimension(3, 3), intent(in) :: eps3 ! 3c epsilon matrix z = r2 - r1, x = dhat
  real(c_double), dimension(3), intent(in) :: r1     ! position of atom 1
  real(c_double), dimension(3), intent(in) :: r2     ! position of atom 2
  real(c_double), dimension(3), intent(in) :: r21    ! vector between atom 1 and atom 2
  real(c_double), dimension(3), intent(in) :: ratm   ! position of the potential
  real(c_double), dimension(3), intent(in) :: rnabc  ! vector from bond charge to potential
  real(c_double), dimension(3, 3, 3), intent(out) :: dera3 ! deps/dratm in the 3-center system
  real(c_double), dimension(3, 3, 3), intent(out) :: der13 ! deps/dr1 in the 3-center system
  integer(c_long) :: imu, ix
  real(c_double) :: crossmag, r1dotr2, r1dotratm, r2dotratm, r2mag2, r1mag2, ratmmag2, sum
  real(c_double), dimension(3) :: crossa
  
  dera3 = 0.0d0
  der13 = 0.0d0
  crossa(1) = r21(2)*rnabc(3) - r21(3)*rnabc(2)
  crossa(2) = r21(3)*rnabc(1) - r21(1)*rnabc(3)
  crossa(3) = r21(1)*rnabc(2) - r21(2)*rnabc(1)
  crossmag = sqrt(crossa(1)*crossa(1) + crossa(2)*crossa(2)  + crossa(3)*crossa(3))
  if (abs(crossmag) .lt. 1.0d-2) then
   return
  end if
  r2mag2 = r2(1)**2 + r2(2)**2 + r2(3)**2
  r1mag2 = r1(1)**2 + r1(2)**2 + r1(3)**2
  ratmmag2 = ratm(1)**2 + ratm(2)**2 + ratm(3)**2
  r1dotr2 = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  r2dotratm = r2(1)*ratm(1) + r2(2)*ratm(2) + r2(3)*ratm(3)
  r1dotratm = r1(1)*ratm(1) + r1(2)*ratm(2) + r1(3)*ratm(3)
  do ix = 1, 3
   do imu = 1, 3
    dera3(ix,imu,3) = 0.0d0
    sum = xlevi(ix,imu,1)*(r2(1) - r1(1))  + xlevi(ix,imu,2)*(r2(2) - r1(2))  + xlevi(ix,imu,3)*(r2(3) - r1(3))
    dera3(ix,imu,2) =  (1.0d0/crossmag)*(sum - (eps3(imu,2)/crossmag)*((r1mag2 + r2mag2 - 2.0d0*r1dotr2)*ratm(ix)  &
    &       + (r1dotr2 - r1mag2 - r2dotratm + r1dotratm)*r2(ix) + (r1dotr2 - r2mag2 + r2dotratm - r1dotratm)*r1(ix)))
   end do
  end do
  dera3(1,1,1) = eps3(3,3)*dera3(1,2,2) - eps3(2,3)*dera3(1,3,2)
  dera3(1,2,1) = eps3(1,3)*dera3(1,3,2) - eps3(3,3)*dera3(1,1,2)
  dera3(1,3,1) = eps3(2,3)*dera3(1,1,2) - eps3(1,3)*dera3(1,2,2)
  dera3(2,1,1) = eps3(3,3)*dera3(2,2,2) - eps3(2,3)*dera3(2,3,2)
  dera3(2,2,1) = eps3(1,3)*dera3(2,3,2) - eps3(3,3)*dera3(2,1,2)
  dera3(2,3,1) = eps3(2,3)*dera3(2,1,2) - eps3(1,3)*dera3(2,2,2)
  dera3(3,1,1) = eps3(3,3)*dera3(3,2,2) - eps3(2,3)*dera3(3,3,2)
  dera3(3,2,1) = eps3(1,3)*dera3(3,3,2) - eps3(3,3)*dera3(3,1,2)
  dera3(3,3,1) = eps3(2,3)*dera3(3,1,2) - eps3(1,3)*dera3(3,2,2)
  do ix = 1, 3
   do imu = 1, 3
    der13(ix,imu,3) = (1.0d0/distance12)*(eps3(imu,3)*eps3(ix,3) - delk(imu,ix))
    sum = xlevi(ix,imu,1)*(ratm(1) - r2(1))  + xlevi(ix,imu,2)*(ratm(2) - r2(2)) + xlevi(ix,imu,3)*(ratm(3) - r2(3))
    der13(ix,imu,2) =  (1.0d0/crossmag)*(sum - (eps3(imu,2)/crossmag)*((r2dotratm - r1dotratm + r1dotr2 - r2mag2)*ratm(ix)       &
    &       + (r2dotratm + r1dotratm - r1dotr2 - ratmmag2)*r2(ix) + (ratmmag2 - 2.0d0*r2dotratm + r2mag2)*r1(ix)))
   end do
  end do
  der13(1,1,1) = eps3(3,3)*der13(1,2,2) - eps3(3,2)*der13(1,2,3) - eps3(2,3)*der13(1,3,2) + eps3(2,2)*der13(1,3,3)
  der13(1,2,1) = eps3(1,3)*der13(1,3,2) - eps3(1,2)*der13(1,3,3) - eps3(3,3)*der13(1,1,2) + eps3(3,2)*der13(1,1,3)
  der13(1,3,1) = eps3(2,3)*der13(1,1,2) - eps3(2,2)*der13(1,1,3) - eps3(1,3)*der13(1,2,2) + eps3(1,2)*der13(1,2,3)
  der13(2,1,1) = eps3(3,3)*der13(2,2,2) - eps3(3,2)*der13(2,2,3) - eps3(2,3)*der13(2,3,2) + eps3(2,2)*der13(2,3,3)
  der13(2,2,1) = eps3(1,3)*der13(2,3,2) - eps3(1,2)*der13(2,3,3) - eps3(3,3)*der13(2,1,2) + eps3(3,2)*der13(2,1,3)
  der13(2,3,1) = eps3(2,3)*der13(2,1,2) - eps3(2,2)*der13(2,1,3) - eps3(1,3)*der13(2,2,2) + eps3(1,2)*der13(2,2,3)
  der13(3,1,1) = eps3(3,3)*der13(3,2,2) - eps3(3,2)*der13(3,2,3) - eps3(2,3)*der13(3,3,2) + eps3(2,2)*der13(3,3,3)
  der13(3,2,1) = eps3(1,3)*der13(3,3,2) - eps3(1,2)*der13(3,3,3) - eps3(3,3)*der13(3,1,2) + eps3(3,2)*der13(3,1,3)
  der13(3,3,1) = eps3(2,3)*der13(3,1,2) - eps3(2,2)*der13(3,1,3) - eps3(1,3)*der13(3,2,2) + eps3(1,2)*der13(3,2,3)
  return
end

