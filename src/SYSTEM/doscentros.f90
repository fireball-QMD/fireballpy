! This subroutine calculates the (two-center) matrix elements (mu,nu).
! There used to be five different routines that did this for all of the
! two-center interactions - doscentros.f, dosxcatm.f, dosxcontop.f,
! dosenatm.f, and dosenontop.f.  These have now all been reduced to one
! routine in order to make Fireball more lean.
!
! This routine also calculates the derivative with respect to the
! position of the orbital of the BRA.
subroutine doscentros (interaction, isub, iforceaux, in1, in2, in3, distance, eps, deps, sx, spx)
  use iso_c_binding
  use M_system
  use M_fdata, only: index_max2c,num_orb,ME2c_max
  implicit none
  integer(c_long), intent (in) :: interaction
  integer(c_long), intent (in) :: isub
  integer(c_long), intent (in) :: iforceaux
  integer(c_long), intent (in) :: in1
  integer(c_long), intent (in) :: in2
  integer(c_long), intent (in) :: in3
  real(c_double), intent (inout) :: distance
  real(c_double), intent (in), dimension (3, 3, 3) :: deps
  real(c_double), intent (in), dimension (3, 3) :: eps
  real(c_double), intent (out), dimension (numorb_max, numorb_max) :: sx
  real(c_double), intent (out), dimension (3, numorb_max, numorb_max) :: spx
  integer(c_long) imu
  integer(c_long) inu
  integer(c_long) index
  real(c_double), dimension (3) :: eta
  real(c_double), dimension (ME2c_max) :: dslist
  real(c_double), dimension (ME2c_max) :: slist
  real(c_double), dimension (numorb_max,numorb_max) :: sm
  real(c_double), dimension (numorb_max,numorb_max) :: spm
  real(c_double), dimension (3,numorb_max,numorb_max) :: spmx
  logical switch

  sm = 0.0d0
  sx = 0.0d0
  if (iforceaux .eq. 1) spm = 0.0d0
  if (iforceaux .eq. 1) spmx = 0.0d0
  if (iforceaux .eq. 1) spx = 0.0d0
  switch = .true.
  if(interaction .eq. 2) switch = .false.
  if(interaction .eq. 15) switch = .false.
  if(interaction .eq. 18) switch = .false.
  do index = 1, index_max2c(in1,in3)
   if ( switch ) then
    call interpolate_1d (interaction, isub, in1, in2, index, iforceaux, distance, slist(index), dslist(index))
   else
    call interpolate_1d (interaction, isub, in1, in3, index, iforceaux, distance, slist(index), dslist(index))
   end if
  end do
  call recover_2c (in1, in3, slist, sm)
  call recover_2c (in1, in3, dslist, spm)
  call rotate_fb (in1, in3, eps, sm, sx)
  if (iforceaux .eq. 1) then
   eta(:) = eps(:,3)
   do inu = 1, num_orb(in3)
    do imu = 1, num_orb(in1)
     if (distance .gt. 1.0d-3) then
      spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
     end if
    end do
   end do
   call rotated (in1, in3, eps, deps, sm, spmx, spx)
  end if
  return
end subroutine doscentros
