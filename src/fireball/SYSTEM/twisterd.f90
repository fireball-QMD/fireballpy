subroutine twisterd (eps, deps, ddmat, dpmat)
  use iso_c_binding
  use M_constants, only: haveDorbitals
  use M_system, only: amat
  implicit none
  real(c_double), dimension(3, 3), intent(in) :: eps
  real(c_double), dimension(3, 3, 3), intent(in) :: deps
  real(c_double), dimension(3, 5, 5), intent(out) :: ddmat
  real(c_double), dimension(3, 3, 3), intent(out) :: dpmat
  integer(c_long) :: imu, ix, jx, kx
  real(c_double) :: aterm12, aterm32, aterm33, aterm13, aterm11, amat_term

  do ix = 1, 3
   dpmat(ix,1,1) = deps(ix,2,2)
   dpmat(ix,1,2) = deps(ix,2,3)
   dpmat(ix,1,3) = deps(ix,2,1)
   dpmat(ix,2,1) = deps(ix,3,2)
   dpmat(ix,2,2) = deps(ix,3,3)
   dpmat(ix,2,3) = deps(ix,3,1)
   dpmat(ix,3,1) = deps(ix,1,2)
   dpmat(ix,3,2) = deps(ix,1,3)
   dpmat(ix,3,3) = deps(ix,1,1)
  end do
  if (.not. haveDorbitals) return
  do imu = 1, 5
   do kx = 1, 3
    aterm12 = 0.0d0
    aterm32 = 0.0d0 
    aterm33 = 0.0d0 
    aterm13 = 0.0d0 
    aterm11 = 0.0d0 
    do ix = 1, 3
     do jx = 1, 3
      if (amat(jx,ix,imu) .ne. 0.0d0) then
       amat_term = amat(jx,ix,imu)
       aterm12 = aterm12 + amat_term*(deps(kx,ix,1)*eps(jx,2) + eps(ix,1)*deps(kx,jx,2))
       aterm32 = aterm32 + amat_term*(deps(kx,ix,3)*eps(jx,2) + eps(ix,3)*deps(kx,jx,2))
       aterm33 = aterm33 + amat_term*(deps(kx,ix,3)*eps(jx,3) + eps(ix,3)*deps(kx,jx,3))
       aterm13 = aterm13 + amat_term*(deps(kx,ix,1)*eps(jx,3) + eps(ix,1)*deps(kx,jx,3))
       aterm11 = aterm11 + amat_term*(deps(kx,ix,1)*eps(jx,1) + eps(ix,1)*deps(kx,jx,1))
      end if
     end do
    end do
    ddmat(kx,imu,1) = 2.0d0*aterm12
    ddmat(kx,imu,2) = 2.0d0*aterm32
    ddmat(kx,imu,3) = sqrt(3.0d0)*aterm33
    ddmat(kx,imu,4) = 2.0d0*aterm13
    ddmat(kx,imu,5) = 2.0d0*aterm11 + aterm33
   end do
  end do
  return
end subroutine twisterd
