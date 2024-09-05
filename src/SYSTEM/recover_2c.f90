subroutine recover_2c (in1, in2, hlist, hbox)
  use iso_c_binding
  use M_system
  use M_fdata, only: ME2c_max, num_orb,index_max2c,mu, num_orb,index_max2c,mu,nu
  implicit none
  integer(c_long), intent(in) :: in1, in2
  real(c_double), intent(in) :: hlist (ME2c_max)
  real(c_double), intent(out) :: hbox (numorb_max, numorb_max)
  integer(c_long) imu, inu
  integer(c_long) index
  do inu = 1, num_orb(in2)
   do imu = 1, num_orb(in1)
    hbox(imu,inu) = 0.0d0
   end do
  end do
  do index = 1, index_max2c(in1,in2)
   imu = mu(index,in1,in2)
   inu = nu(index,in1,in2)
   hbox(imu,inu) = hlist (index)
  end do
  return
end

