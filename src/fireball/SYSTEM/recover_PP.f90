subroutine recover_PP (in1, in2, hlist, hbox)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max
  use M_fdata, only: ME2c_max,num_orbPP,num_orb,index_maxPP,nuPP,muPP
  implicit none
  integer, intent(in) :: in1, in2
  real(double), intent(in), dimension (ME2c_max) :: hlist
  real(double), intent(out), dimension (numorb_max, numorb_max) :: hbox
  integer imu, inu
  integer index
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
