subroutine rotatePP(in1,in2,eps,mmatrix,xmatrix)
  use iso_c_binding
  use M_system
  use M_fdata, only: nssh,lssh,nsshPP,lsshPP
  implicit none
  integer(c_long), intent(in) :: in1
  integer(c_long), intent(in) :: in2
  real(c_double), intent(in) :: eps (3, 3)
  real(c_double), intent(in) :: mmatrix (numorb_max, numorb_max)
  real(c_double), intent(out) :: xmatrix (numorb_max, numorb_max)
  integer(c_long) issh
  integer(c_long) jssh
  integer(c_long) k1, k2
  integer(c_long) n1, l1, m1
  integer(c_long) n2, l2, m2
  real(c_double) dmat (5, 5)
  real(c_double) left (5, 5)
  real(c_double) pmat (3, 3)
  real(c_double) right (5, 5)
  call twister (eps, dmat, pmat)
  xmatrix=0.0d0
  n1 = 0
  do issh = 1, nssh(in1)
   l1 = lssh(issh,in1)
   call chooser (l1, dmat, pmat, left)
   n2 = 0
   do jssh = 1, nsshPP(in2)
    l2 = lsshPP(jssh,in2)
    call chooser (l2, dmat, pmat, right)
    do m2 = 1, 2*l2 + 1
     do m1 = 1, 2*l1 + 1
      do k2 = 1, 2*l2 + 1
       do k1 = 1, 2*l1 + 1
        xmatrix(n1+k1,n2+k2) = xmatrix(n1+k1,n2+k2)   &
     &    + left(k1,m1)*mmatrix(n1+m1,n2+m2)*right(k2,m2)
       end do
      end do
     end do
    end do
    n2 = n2 + 2*l2 + 1
   end do
   n1 = n1 + 2*l1 + 1
  end do
  return
end subroutine rotatePP
