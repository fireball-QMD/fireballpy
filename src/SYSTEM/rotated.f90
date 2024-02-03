subroutine rotated (in1, in2, eps, deps, matm, dmatm, dmatx)
  use M_system
  implicit none
  integer, intent(in) :: in1, in2
  real, intent(in) :: deps (3, 3, 3)
  real, intent(in) :: eps (3, 3)
  real, intent(in) :: dmatm (3, numorb_max, numorb_max)
  real, intent(in) :: matm (numorb_max, numorb_max)
  real, intent(out) :: dmatx (3, numorb_max, numorb_max)
  integer imu
  integer inu
  integer ix
  real ddmat (3, 5, 5)
  real dmat (5, 5)
  real dpmat (3, 3, 3)
  real pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end

