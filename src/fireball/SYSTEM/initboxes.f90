  subroutine initboxes ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: mbeta_max, xl, a1vec, a2vec, a3vec
  implicit none
  integer lbeta
  integer ix, iy
  integer mbeta
  integer midl
  real(dp), save, dimension (3, 0:728) :: xxl ! cube 9x9x9
  integer, parameter :: mbox = 4
  xxl (:,:) = 0.0_dp
  xxl (1,1) = 1.0_dp
  xxl (1,2) = -1.0_dp
  xxl (2,3) = 1.0_dp
  xxl (2,4) = -1.0_dp
  xxl (3,5) = 1.0_dp
  xxl (3,6) = -1.0_dp
  xxl (1:2,7) = 1.0_dp
  xxl (1,8) = -1.0_dp
  xxl (2,8) = 1.0_dp
  xxl (1:2,9) = -1.0_dp
  xxl (1,10) = 1.0_dp
  xxl (2,10) = -1.0_dp
  xxl (:,11) = 1.0_dp
  xxl (1,12) = -1.0_dp
  xxl (2:3,12) = 1.0_dp
  xxl (1:2,13) = -1.0_dp
  xxl (3,13) = 1.0_dp
  xxl (:,14) = 1.0_dp
  xxl (2,14) = -1.0_dp
  xxl (1:2,15) = 1.0_dp
  xxl (3,15) = -1.0_dp
  xxl (:,16) = -1.0_dp
  xxl (2,16) = 1.0_dp
  xxl (:,17) = -1.0_dp
  xxl (1,18) = 1.0_dp
  xxl (2:3,18) = -1.0_dp
  xxl (1,19) = 1.0_dp
  xxl (3,19) = 1.0_dp
  xxl (1,20) = -1.0_dp
  xxl (3,20) = 1.0_dp
  xxl (2:3,21) = 1.0_dp
  xxl (2,22) = -1.0_dp
  xxl (3,22) = 1.0_dp
  xxl (1,23) = 1.0_dp
  xxl (3,23) = -1.0_dp
  xxl (1,24) = -1.0_dp
  xxl (3,24) = -1.0_dp
  xxl (2,25) = 1.0_dp
  xxl (3,25) = -1.0_dp
  xxl (2:3,26) = -1.0_dp
  lbeta = 26
  do ix = -mbox, mbox
    do iy = -mbox, mbox
      lbeta = lbeta + 1
      xxl(1,lbeta) = real(ix, dp)
      xxl(2,lbeta) = real(iy, dp)
      xxl(3,lbeta) = -4.0_dp
    end do
  end do
  do ix = -mbox, mbox
    do iy = -mbox, mbox
      lbeta = lbeta + 1
      xxl(1,lbeta) = real(ix, dp)
      xxl(2,lbeta) = real(iy, dp)
      xxl(3,lbeta) = -3.0_dp
    end do
  end do
  do ix = -mbox, mbox
    do iy = -mbox, mbox
      lbeta = lbeta + 1
      xxl(1,lbeta) = real(ix, dp)
      xxl(2,lbeta) = real(iy, dp)
      xxl(3,lbeta) = -2.0_dp
    end do
  end do
  do midl = -1, 1
   do ix = -mbox, mbox
     do iy = -mbox, mbox
       if (abs(ix) .gt. 1 .or. abs(iy) .gt. 1) then
        lbeta = lbeta + 1
        xxl(1,lbeta) = real(ix, dp)
        xxl(2,lbeta) = real(iy, dp)
        xxl(3,lbeta) = real(midl, dp)
       end if
     end do
   end do
  end do
  do ix = -mbox, mbox
   do iy = -mbox, mbox
     lbeta = lbeta + 1
     xxl(1,lbeta) = real(ix, dp)
     xxl(2,lbeta) = real(iy, dp)
     xxl(3,lbeta) = 2.0_dp
   end do
  end do
  do ix = -mbox, mbox
  do iy = -mbox, mbox
    lbeta = lbeta + 1
     xxl(1,lbeta) = real(ix, dp)
     xxl(2,lbeta) = real(iy, dp)
     xxl(3,lbeta) = 3.0_dp
   end do
  end do
  do ix = -mbox, mbox
    do iy = -mbox, mbox
      lbeta = lbeta + 1
      xxl(1,lbeta) = real(ix, dp)
      xxl(2,lbeta) = real(iy, dp)
       xxl(3,lbeta) = 4.0_dp
     end do
    end do
  mbeta_max = lbeta
  if (allocated(xl)) deallocate(xl)
  allocate ( xl(3,0:mbeta_max) )
  do mbeta = 0,mbeta_max
    xl(:,mbeta) = xxl(1,mbeta)*a1vec(:) + xxl(2,mbeta)*a2vec(:)  + xxl(3,mbeta)*a3vec(:)
  end do
  return
end subroutine initboxes
