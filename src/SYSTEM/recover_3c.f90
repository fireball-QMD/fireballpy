subroutine recover_3c (in1, in2, hlist, hbox)
  use iso_c_binding
  use M_system, only: numorb_max
  use m_fdata, only: ME3c_max,num_orb,index_max3c,mu,nu
  implicit none
  integer(c_long), intent(in) :: in1, in2
  real(c_double), intent(in) :: hlist (ME3c_max)
  real(c_double), intent(out) :: hbox (numorb_max, numorb_max)
  integer(c_long) imu, inu
  integer(c_long) index
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
end subroutine recover_3c
