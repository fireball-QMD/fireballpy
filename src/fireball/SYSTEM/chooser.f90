subroutine chooser (l, dmat, pmat, rmatrix)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: l
  real(dp), intent(in), dimension(5, 5) :: dmat
  real(dp), intent(in), dimension(3, 3) :: pmat
  real(dp), intent(out), dimension(5, 5) :: rmatrix
  rmatrix = 0.0_dp
  if (l .eq. 0) then
    rmatrix(1,1) = 1.0_dp
  else if (l .eq. 1) then
    rmatrix(1:3,1:3) = pmat
  else if (l .eq. 2) then
    rmatrix = dmat
  end if
  return
end subroutine chooser
