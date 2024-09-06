subroutine doscentrosDipY (interaction, isub, in1, in2, in3, distance, eps, deps, sx, spx)
  use iso_c_binding
  use M_system
  use M_fdata, only: num_orb, index_max2cDipY,ME2cDipY_max
  implicit none
  integer(c_long), intent (in) :: interaction
  integer(c_long), intent (in) :: isub
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
  real(c_double), dimension (ME2cDipY_max) :: dslist
  real(c_double), dimension (ME2cDipY_max) :: slist
  real(c_double), dimension (numorb_max,numorb_max) :: sm
  real(c_double), dimension (numorb_max,numorb_max) :: spm
  real(c_double), dimension (3,numorb_max,numorb_max) :: spmx
  sm = 0.0d0
  sx = 0.0d0
  if (iforce .eq. 1) spm = 0.0d0
  if (iforce .eq. 1) spmx = 0.0d0
  if (iforce .eq. 1) spx = 0.0d0
  do index = 1, index_max2cDipY(in1,in3)
    call interpolate_1d (interaction, isub, in1, in2, index, iforce, distance, slist(index), dslist(index))
  end do
  call recover_2cDipY (in1, in3, slist, sm)
  call recover_2cDipY (in1, in3, dslist, spm)
  call rotate_fb (in1, in3, eps, sm, sx)
  if (iforce .eq. 1) then
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
end subroutine doscentrosDipY
