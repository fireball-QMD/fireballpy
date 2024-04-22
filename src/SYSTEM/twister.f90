subroutine twister (eps, dmat, pmat)
  use M_constants
  implicit none
  real(8), intent(in) :: eps (3, 3)
  real(8), intent(out) :: dmat (5, 5)
  real(8), intent(out) :: pmat (3, 3)
  integer imu
  integer jx
  integer ix
  real(8) amat_term
  real(8) xlambda11
  real(8) xlambda12
  real(8) xlambda13
  real(8) xlambda32
  real(8) xlambda33
  pmat(1,1) = eps(2,2) 
  pmat(1,2) = eps(2,3) 
  pmat(1,3) = eps(2,1)
  pmat(2,1) = eps(3,2) 
  pmat(2,2) = eps(3,3)
  pmat(2,3) = eps(3,1)
  pmat(3,1) = eps(1,2)
  pmat(3,2) = eps(1,3)
  pmat(3,3) = eps(1,1)
  if (.not. haveDorbitals) return
  do imu = 1, 5
   xlambda11 = 0
   xlambda12 = 0
   xlambda13 = 0
   xlambda32 = 0
   xlambda33 = 0
   do jx = 1, 3
    do ix = 1, 3
     if (amat(ix,jx,imu) .ne. 0.0d0) then
      amat_term = amat(ix,jx,imu)
      xlambda11 = xlambda11 + amat_term*eps(jx,1)*eps(ix,1)
      xlambda12 = xlambda12 + amat_term*eps(jx,1)*eps(ix,2)
      xlambda13 = xlambda13 + amat_term*eps(jx,1)*eps(ix,3)
      xlambda32 = xlambda32 + amat_term*eps(jx,3)*eps(ix,2)
      xlambda33 = xlambda33 + amat_term*eps(jx,3)*eps(ix,3)
     end if
    end do
   end do
   dmat(imu,1) = 2.0d0*xlambda12
   dmat(imu,2) = 2.0d0*xlambda32
   dmat(imu,3) = sqrt(3.0d0)*xlambda33
   dmat(imu,4) = 2.0d0*xlambda13
   dmat(imu,5) = 2.0d0*xlambda11 + xlambda33
  end do
  return
end

