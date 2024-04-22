subroutine rotated (in1, in2, eps, deps, matm, dmatm, dmatx)
  use M_system
  implicit none
  integer, intent(in) :: in1, in2
  real(8), intent(in) :: deps (3, 3, 3)
  real(8), intent(in) :: eps (3, 3)
  real(8), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(8), intent(in) :: matm (numorb_max, numorb_max)
  real(8), intent(out) :: dmatx (3, numorb_max, numorb_max)
  integer imu
  integer inu
  integer ix
  real(8) ddmat (3, 5, 5)
  real(8) dmat (5, 5)
  real(8) dpmat (3, 3, 3)
  real(8) pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end

