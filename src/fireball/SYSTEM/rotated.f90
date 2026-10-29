subroutine rotated (in1, in2, eps, deps, matm, dmatm, dmatx)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: numorb_max
  implicit none
  integer, intent(in) :: in1, in2
  real(dp), intent(in) :: deps (3, 3, 3)
  real(dp), intent(in) :: eps (3, 3)
  real(dp), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(dp), intent(in) :: matm (numorb_max, numorb_max)
  real(dp), intent(out) :: dmatx (3, numorb_max, numorb_max)
  real(dp) ddmat (3, 5, 5)
  real(dp) dmat (5, 5)
  real(dp) dpmat (3, 3, 3)
  real(dp) pmat (3, 3)
  call twister (eps, dmat, pmat)
  call twisterd (eps, deps, ddmat, dpmat)
  call makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, dmatx)
  return
end subroutine rotated
