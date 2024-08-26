subroutine recover_S (in1, in2, hlist, hbox)
  use M_constants, only: wp
  use M_system
  use M_fdata, only: nssh, index_maxS, muS, nuS, MES_max, nsh_max, nsh_max
  implicit none
  integer, intent(in) :: in1, in2
  real(wp), intent(in) :: hlist (MES_max)
  real(wp), intent(out) :: hbox (nsh_max, nsh_max)
  integer imu, inu
  integer index
  do inu = 1, nssh(in2)
    do imu = 1, nssh(in1)
      hbox(imu,inu) = 0.0d0
    end do
  end do
  do index = 1, index_maxS(in1,in2)
    imu = muS(index,in1,in2)
    inu = nuS(index,in1,in2)
    hbox(imu,inu) = hlist (index)
  end do
  return
end subroutine recover_S

