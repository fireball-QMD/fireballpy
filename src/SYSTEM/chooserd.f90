subroutine chooserd (l, ddmat, dpmat, dmatrix)
  implicit none
  integer, intent (in) :: l
  real, intent (in), dimension(3, 5, 5) :: ddmat
  real, intent (in), dimension(3, 3, 3) :: dpmat
  real, intent (out), dimension(3, 5, 5) :: dmatrix
  dmatrix = 0.0d0
  if (l .eq. 0) then
  else if (l .eq. 1) then
   dmatrix(:,1:3,1:3) = dpmat
  else if (l .eq. 2) then
   dmatrix = ddmat
  end if
  return
end

