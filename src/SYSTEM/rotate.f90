subroutine rotate_fb (in1, in2, eps, mmatrix, xmatrix)
  use M_constants, only: wp
  use M_system
  use M_fdata, only: nssh,lssh
  implicit none
  integer, intent(in) :: in1
  integer, intent(in) :: in2
  real(wp), intent(in) :: eps (3, 3)
  real(wp), intent(in) :: mmatrix (numorb_max, numorb_max)
  real(wp), intent(out) :: xmatrix (numorb_max, numorb_max)
  integer issh
  integer jssh
  integer k1, k2
  integer n1, l1, m1
  integer n2, l2, m2
  real(wp) dmat (5, 5)
  real(wp) left (5, 5)
  real(wp) pmat (3, 3)
  real(wp) right (5, 5)
  call twister (eps, dmat, pmat)
  xmatrix=0.0d0
  n1 = 0
  do issh = 1, nssh(in1)
   l1 = lssh(issh,in1)
   call chooser (l1, dmat, pmat, left)
   n2 = 0
   do jssh = 1, nssh(in2)
    l2 = lssh(jssh,in2)
    call chooser (l2, dmat, pmat, right)
    do m2 = 1, 2*l2 + 1
     do m1 = 1, 2*l1 + 1
      do k2 = 1, 2*l2 + 1
       do k1 = 1, 2*l1 + 1
        xmatrix(n1+k1,n2+k2) = xmatrix(n1+k1,n2+k2) + left(k1,m1)*mmatrix(n1+m1,n2+m2)*right(k2,m2)
       end do
      end do
     end do
    end do
    n2 = n2 + 2*l2 + 1
   end do
   n1 = n1 + 2*l1 + 1
  end do
  return
end

