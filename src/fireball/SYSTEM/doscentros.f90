! This subroutine calculates the (two-center) matrix elements (mu,nu).
! There used to be five different routines that did this for all of the
! two-center interactions - doscentros.f, dosxcatm.f, dosxcontop.f,
! dosenatm.f, and dosenontop.f.  These have now all been reduced to one
! routine in order to make Fireball more lean.
!
! This routine also calculates the derivative with respect to the
! position of the orbital of the BRA.
subroutine doscentros (interaction, isub, iforceaux, in1, in2, in3, distance, eps, deps, sx, spx)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: numorb_max
  use M_fdata, only: index_max2c,num_orb,ME2c_max, TWOCENTER_VNA_L, TWOCENTER_VXC_L, TWOCENTER_DEN_L, splineint_2c, ind2c
  implicit none
  integer, intent (in) :: interaction
  integer, intent (in) :: isub
  integer, intent (in) :: iforceaux
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
  real(dp), dimension (ME2c_max) :: dslist
  real(dp), dimension (ME2c_max) :: slist
  real(dp), dimension (numorb_max,numorb_max) :: sm
  real(dp), dimension (numorb_max,numorb_max) :: spm
  real(dp), dimension (3,numorb_max,numorb_max) :: spmx
  logical switch

  sm = 0.0_dp
  sx = 0.0_dp
  if (iforceaux .eq. 1) spm = 0.0_dp
  if (iforceaux .eq. 1) spmx = 0.0_dp
  if (iforceaux .eq. 1) spx = 0.0_dp
  switch = .true.
  if(interaction .eq. TWOCENTER_VNA_L) switch = .false.
  if(interaction .eq. TWOCENTER_VXC_L) switch = .false.
  if(interaction .eq. TWOCENTER_DEN_L) switch = .false.
  do index = 1, index_max2c(in1,in3)
   if ( switch ) then
    slist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance)
    if (iforceaux == 1) dslist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance, order=1)
   else
    slist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in3)%f(distance)
    if (iforceaux == 1) dslist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in3)%f(distance, order=1)
   end if
  end do
  call recover_2c (in1, in3, slist, sm)
  call recover_2c (in1, in3, dslist, spm)
  call rotate_fb (in1, in3, eps, sm, sx)
  if (iforceaux .eq. 1) then
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
end subroutine doscentros
