! This routine computes the non-local pseudopotential matrix element ppx(mu,nu)
! The range of neighbours could/should be different for  3-center pseudopotential matrix elements than for neutral-coulomb contributions.  This could/should be taken into account in the future
subroutine doscentrosPP (interaction, isub, distance, eps, deps, iauxforce, in1, in2, sx, spx)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: numorb_max
  use M_fdata, only: index_maxPP,num_orb,num_orbpp,ME2cPP_max,splineint_2c,ind2c
  implicit none
  integer, intent(in) :: iauxforce
  integer, intent(in) :: in1
  integer, intent(in) :: in2
  integer, intent(in) :: isub
  integer, intent(in) :: interaction
  real(dp), intent (in) :: distance
  real(dp), intent (in), dimension (3, 3, 3) :: deps
  real(dp), intent (in), dimension (3, 3) :: eps
  real(dp), intent(out), dimension (numorb_max, numorb_max) :: sx
  real(dp), intent(out), dimension (3, numorb_max, numorb_max) :: spx
  integer imu
  integer inu
  integer index
  real(dp), dimension (ME2cPP_max) :: dpplist
  real(dp), dimension (3) :: eta
  real(dp), dimension (ME2cPP_max) :: pplist
  real(dp), dimension (numorb_max, numorb_max) :: sm
  real(dp), dimension (numorb_max, numorb_max) :: spm
  real(dp), dimension (3, numorb_max, numorb_max) :: spmx
  sm = 0.0_dp
  sx = 0.0_dp
  if (iauxforce .eq. 1) spx = 0.0_dp
  do index = 1, index_maxPP(in1,in2)
    pplist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance)
    if (iauxforce == 1) dpplist(index) = splineint_2c(index, ind2c(interaction, isub), in1, in2)%f(distance, order=1)
  end do
  call recover_PP (in1, in2, pplist, sm)
  call recover_PP (in1, in2, dpplist, spm)
  call rotatePP (in1, in2, eps, sm, sx)
  if (iauxforce .eq. 1) then
    eta(:) = eps(:,3)
    do imu = 1, num_orb(in1)
      do inu = 1, num_orbpp(in2)
        spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
      end do
    end do
    call rotatedPP (in1, in2, eps, deps, sm, spmx, spx)
  end if
  return
end subroutine doscentrosPP
