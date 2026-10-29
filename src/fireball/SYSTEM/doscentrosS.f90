subroutine doscentrosS (interaction, isub, iauxforce, in1, in2, in3, distance, eps, sx, spx)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_fdata, only: index_maxS,nsh_max,MES_max,nssh, TWOCENTER_DENS_L,splineint_2c,ind2c
  implicit none
  integer, intent (in) :: iauxforce
  integer, intent (in) :: interaction
  integer, intent (in) :: isub
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: in3
  real(dp), intent (inout) :: distance
  real(dp), intent (in), dimension (3, 3) :: eps
  real(dp), intent (out), dimension (nsh_max, nsh_max) :: sx
  real(dp), intent (out), dimension (3, nsh_max, nsh_max) :: spx
  integer imu
  integer inu
  integer index
  real(dp), dimension (3) :: eta
  real(dp), dimension (MES_max) :: dslist
  real(dp), dimension (MES_max) :: slist
  real(dp), dimension (nsh_max,nsh_max) :: spm

  sx = 0.0_dp
  if (iauxforce .eq. 1) spm = 0.0_dp
  if (iauxforce .eq. 1) spx = 0.0_dp
  do index = 1, index_maxS(in1,in3)
    if (interaction .ne. TWOCENTER_DENS_L) then
      slist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance)
      if (iauxforce == 1) dslist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance, order=1)
    else
      slist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in3)%f(distance)
      if (iauxforce == 1) dslist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in3)%f(distance, order=1)
    end if
  end do
  call recover_S (in1, in3, slist, sx)
  call recover_S (in1, in3, dslist, spm)
  if (iauxforce .eq. 1) then
    eta(:) = eps(:,3)
    do inu = 1, nssh(in3)
      do imu = 1, nssh(in1)
        if (distance .gt. 1.0e-3_dp) then
          spx(:,imu,inu) = - eta(:)*spm(imu,inu)
        end if
      end do
    end do
  end if
  return
end subroutine doscentrosS
