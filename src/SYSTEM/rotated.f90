subroutine rotated (in1, in2, eps, deps, matm, dmatm, dmatx)
  use iso_c_binding
  use M_system, only: numorb_max
  implicit none
  integer(c_long), intent(in) :: in1, in2
  real(c_double), intent(in) :: deps (3, 3, 3)
  real(c_double), intent(in) :: eps (3, 3)
  real(c_double), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(c_double), intent(in) :: matm (numorb_max, numorb_max)
  real(c_double), intent(out) :: dmatx (3, numorb_max, numorb_max)
  real(c_double) ddmat (3, 5, 5)
  real(c_double) dmat (5, 5)
  real(c_double) dpmat (3, 3, 3)
  real(c_double) pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end subroutine rotated
