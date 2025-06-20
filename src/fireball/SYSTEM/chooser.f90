subroutine chooser (l, dmat, pmat, rmatrix)
  use, intrinsic :: iso_fortran_env, only: double => real64
  implicit none
  integer, intent(in) :: l
  real(double), intent(in), dimension(5, 5) :: dmat
  real(double), intent(in), dimension(3, 3) :: pmat
  real(double), intent(out), dimension(5, 5) :: rmatrix
  rmatrix = 0.0d0
  if (l .eq. 0) then
    rmatrix(1,1) = 1.0d0
  else if (l .eq. 1) then
    rmatrix(1:3,1:3) = pmat
  else if (l .eq. 2) then
    rmatrix = dmat
  end if
  return
end subroutine chooser
