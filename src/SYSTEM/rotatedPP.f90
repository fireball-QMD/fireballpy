subroutine rotatedPP (in1, in2, eps, deps, matm, dmatm, dmatx)
  use M_constants, only: wp
  use M_system
  implicit none
  integer, intent(in) :: in1, in2
  real(wp), intent(in) :: deps (3, 3, 3)
  real(wp), intent(in) :: eps (3, 3)
  real(wp), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(wp), intent(in) :: matm (numorb_max, numorb_max)
  real(wp), intent(out) :: dmatx (3, numorb_max, numorb_max)
  integer imu
  integer inu
  integer ix
  real(wp) ddmat (3, 5, 5)
  real(wp) dmat (5, 5)
  real(wp) dpmat (3, 3, 3)
  real(wp) pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmatPP (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end

