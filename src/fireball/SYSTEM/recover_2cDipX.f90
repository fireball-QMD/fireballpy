subroutine recover_2cDipX (in1, in2, hlist, hbox)
  use iso_c_binding
  use M_system, only: numorb_max
  use M_fdata, only: num_orb, muDipX, nuDipX, index_max2cDipX,ME2cDipX_max
  implicit none
  integer(c_long), intent(in) :: in1, in2
  real(c_double), intent(in) :: hlist (ME2cDipX_max)
  real(c_double), intent(out) :: hbox (numorb_max, numorb_max)
  integer(c_long) imu, inu
  integer(c_long) index
  do inu = 1, num_orb(in2)
   do imu = 1, num_orb(in1)
    hbox(imu,inu) = 0.0d0
   end do
  end do
  do index = 1, index_max2cDipX(in1,in2)
   imu = muDipX(index,in1,in2)
   inu = nuDipX(index,in1,in2)
   hbox(imu,inu) = hlist (index)
  end do
  return
end subroutine recover_2cDipX
