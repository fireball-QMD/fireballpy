! This routine computes the non-local pseudopotential matrix element ppx(mu,nu)
! The range of neighbours could/should be different for  3-center pseudopotential matrix elements than for neutral-coulomb contributions.  This could/should be taken into account in the future
subroutine doscentrosPP (interaction, isub, distance, eps, deps, iauxforce, in1, in2, sx, spx)
  use iso_c_binding
  use M_system
  use M_fdata, only: index_maxPP,num_orb,num_orbpp,ME2cPP_max
  implicit none
  integer(c_long), intent(in) :: iauxforce
  integer(c_long), intent(in) :: in1
  integer(c_long), intent(in) :: in2
  integer(c_long), intent(in) :: isub
  integer(c_long), intent(in) :: interaction
  real(c_double), intent (in) :: distance
  real(c_double), intent (in), dimension (3, 3, 3) :: deps
  real(c_double), intent (in), dimension (3, 3) :: eps
  real(c_double), intent(out), dimension (numorb_max, numorb_max) :: sx
  real(c_double), intent(out), dimension (3, numorb_max, numorb_max) :: spx
  integer(c_long) imu
  integer(c_long) inu
  integer(c_long) index
  real(c_double), dimension (ME2cPP_max) :: dpplist
  real(c_double), dimension (3) :: eta
  real(c_double), dimension (ME2cPP_max) :: pplist
  real(c_double), dimension (numorb_max, numorb_max) :: sm
  real(c_double), dimension (numorb_max, numorb_max) :: spm
  real(c_double), dimension (3, numorb_max, numorb_max) :: spmx
  if (interaction .ne. 5 .or. isub .ne. 0) then
    write (*,*) ' interaction = ', interaction
    write (*,*) ' This routine is only for the pseudopotential '
    write (*,*) ' interactions.  The wrong subroutine is being called. '
    write (*,*) ' We must stop here! '
    stop
  end if
  sm = 0.0d0
  sx = 0.0d0
  if (iauxforce .eq. 1) spx = 0.0d0
  do index = 1, index_maxPP(in1,in2)
    call interpolate_1d (interaction, isub, in1, in2, index, iauxforce, distance, pplist(index), dpplist(index))
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
end
