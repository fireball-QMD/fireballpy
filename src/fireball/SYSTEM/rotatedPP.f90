subroutine rotatedPP (in1, in2, eps, deps, matm, dmatm, dmatx)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max
  implicit none
  integer, intent(in) :: in1, in2
  real(double), intent(in) :: deps (3, 3, 3)
  real(double), intent(in) :: eps (3, 3)
  real(double), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(double), intent(in) :: matm (numorb_max, numorb_max)
  real(double), intent(out) :: dmatx (3, numorb_max, numorb_max)
  real(double) ddmat (3, 5, 5)
  real(double) dmat (5, 5)
  real(double) dpmat (3, 3, 3)
  real(double) pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmatPP (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end subroutine rotatedPP
