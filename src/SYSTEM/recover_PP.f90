subroutine recover_PP (in1, in2, hlist, hbox)
  use iso_c_binding
  use M_system, only: numorb_max
  use M_fdata, only: ME2c_max,num_orbPP,num_orb,index_maxPP,nuPP,muPP
  implicit none
  integer(c_long), intent(in) :: in1, in2
  real(c_double), intent(in), dimension (ME2c_max) :: hlist
  real(c_double), intent(out), dimension (numorb_max, numorb_max) :: hbox
  integer(c_long) imu, inu
  integer(c_long) index
  do inu = 1, num_orbPP(in2)
   do imu = 1, num_orb(in1)
    hbox(imu,inu) = 0.0d0
   end do
  end do
  do index = 1, index_maxPP(in1,in2)
   imu = muPP(index,in1,in2)
   inu = nuPP(index,in1,in2)
   hbox(imu,inu) = hlist(index)
  end do
  return
end subroutine recover_PP
