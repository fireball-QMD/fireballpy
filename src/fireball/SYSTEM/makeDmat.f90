subroutine makeDmat (in1, in2, matm, dmatm, dmat, pmat, ddmat, dpmat, term)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max
  use M_fdata, only: nssh,lssh
  implicit none
  integer, intent(in) :: in1, in2
  real(double), intent(in) :: dmat (5, 5)
  real(double), intent(in) :: dmatm (3, numorb_max, numorb_max)
  real(double), intent(in) :: ddmat (3, 5, 5)
  real(double), intent(in) :: matm (numorb_max, numorb_max)
  real(double), intent(in) :: pmat (3, 3)
  real(double), intent(in) :: dpmat (3, 3, 3)
  real(double), intent(out) :: term (3, numorb_max, numorb_max)
  integer issh, jssh
  integer ix
  integer k1, k2
  integer l1, l2
  integer m1, m2
  integer n1, n2
  real(double) dleft (3, 5, 5)
  real(double) dright (3, 5, 5)
  real(double) left (5, 5)
  real(double) right (5, 5)
  n1 = 0
  do issh = 1, nssh(in1)
   l1 = lssh(issh,in1)
   call chooserd (l1, ddmat, dpmat, dleft)
   call chooser (l1, dmat, pmat, left)
   n2 = 0
   do jssh = 1, nssh(in2)
    l2 = lssh(jssh,in2)
    call chooserd (l2, ddmat, dpmat, dright)
    call chooser (l2, dmat, pmat, right)
    do k2 = 1, 2*l2 + 1
     do k1 = 1, 2*l1 + 1
      do ix = 1, 3
       term(ix,n1+k1,n2+k2) = 0.0d0
      end do
     end do
    end do
    do m2 = 1, 2*l2 + 1
     do k2 = 1, 2*l2 + 1
      do m1 = 1, 2*l1 + 1
       do k1 = 1, 2*l1 + 1
        do ix = 1, 3
          term(ix,n1+k1,n2+k2) = term(ix,n1+k1,n2+k2)  &
          &     + dleft(ix,k1,m1)*right(k2,m2)    *matm(n1+m1,n2+m2)    &
          &     + left(k1,m1)    *dright(ix,k2,m2)*matm(n1+m1,n2+m2)    &
          &     + left(k1,m1)    *right(k2,m2)    *dmatm(ix,n1+m1,n2+m2)
        end do
       end do
      end do
     end do
    end do
    n2 = n2 + 2*l2 + 1
   end do
   n1 = n1 + 2*l1 + 1
  end do
  return
end subroutine makeDmat
