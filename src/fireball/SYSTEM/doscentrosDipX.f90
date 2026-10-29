subroutine doscentrosDipX (interaction, isub, in1, in2, in3, distance, eps, deps, sx, spx)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: iforce, numorb_max
  use M_fdata, only: num_orb,index_max2cDipX,ME2cDipX_max,splineint_2c,ind2c
  implicit none
  integer, intent (in) :: interaction
  integer, intent (in) :: isub
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: in3
  real(dp), intent (inout) :: distance
  real(dp), intent (in), dimension (3, 3, 3) :: deps
  real(dp), intent (in), dimension (3, 3) :: eps
  real(dp), intent (out), dimension (numorb_max, numorb_max) :: sx
  real(dp), intent (out), dimension (3, numorb_max, numorb_max) :: spx
  integer imu
  integer inu
  integer index
  real(dp), dimension (3) :: eta
  real(dp), dimension (ME2cDipX_max) :: dslist
  real(dp), dimension (ME2cDipX_max) :: slist
  real(dp), dimension (numorb_max,numorb_max) :: sm
  real(dp), dimension (numorb_max,numorb_max) :: spm
  real(dp), dimension (3,numorb_max,numorb_max) :: spmx
  sm = 0.0_dp
  sx = 0.0_dp
  if (iforce .eq. 1) spm = 0.0_dp
  if (iforce .eq. 1) spmx = 0.0_dp
  if (iforce .eq. 1) spx = 0.0_dp
  do index = 1, index_max2cDipX(in1,in3)
    slist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance)
    if (iforce == 1) dslist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance, order=1)
  end do
  call recover_2cDipX (in1, in3, slist, sm)
  call recover_2cDipX (in1, in3, dslist, spm)
  call rotate_fb (in1, in3, eps, sm, sx)
  if (iforce .eq. 1) then
    eta(:) = eps(:,3)
    do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
        if (distance .gt. 1.0e-3_dp) then
          spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
        end if
      end do
    end do
    call rotated (in1, in3, eps, deps, sm, spmx, spx)
  end if
  return
end subroutine doscentrosDipX
