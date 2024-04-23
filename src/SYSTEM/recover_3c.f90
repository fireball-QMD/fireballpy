subroutine recover_3c (in1, in2, hlist, hbox)
  use M_system
  use m_fdata, only: ME3c_max,num_orb,index_max3c,mu,nu
  implicit none
  integer, intent(in) :: in1, in2
  real*8, intent(in) :: hlist (ME3c_max)
  real*8, intent(out) :: hbox (numorb_max, numorb_max)
  integer imu, inu
  integer index
  do inu = 1, num_orb(in2)
   do imu = 1, num_orb(in1)
    hbox(imu,inu) = 0.0d0
   end do
  end do
  do index = 1, index_max3c(in1,in2)
   imu = mu(index,in1,in2)
   inu = nu(index,in1,in2)
   hbox(imu,inu) = hlist(index)
  end do
  return
end

