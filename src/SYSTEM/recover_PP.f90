! This subroutine calculates the BOX Hbox (num_orb(in1) x num_orb(in2) of matrix-elements from the list of matrix elements stored in Hlist.
subroutine recover_PP (in1, in2, hlist, hbox)
  use M_system
  implicit none
  integer, intent(in) :: in1, in2
  real, intent(in), dimension (ME2c_max) :: hlist
  real, intent(out), dimension (numorb_max, numorb_max) :: hbox
 
  integer imu, inu
  integer index
 
! Initialize hbox
  do inu = 1, num_orbPP(in2)
   do imu = 1, num_orb(in1)
    hbox(imu,inu) = 0.0d0
   end do
  end do
 
! Now, construct hbox
  do index = 1, index_maxPP(in1,in2)
   imu = muPP(index,in1,in2)
   inu = nuPP(index,in1,in2)
   hbox(imu,inu) = hlist(index)
  end do
 
  return
end
