subroutine doscentrosDipY (interaction, isub, in1, in2, in3, distance, eps, deps, sx, spx)
  use M_system
  use M_fdata, only: num_orb, index_max2cDipY,ME2cDipY_max
  implicit none
  integer, intent (in) :: interaction
  integer, intent (in) :: isub
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: in3
  real(8), intent (inout) :: distance
  real(8), intent (in), dimension (3, 3, 3) :: deps
  real(8), intent (in), dimension (3, 3) :: eps
  real(8), intent (out), dimension (numorb_max, numorb_max) :: sx
  real(8), intent (out), dimension (3, numorb_max, numorb_max) :: spx

  integer imu
  integer inu
  integer index
  real(8), dimension (3) :: eta
  real(8), dimension (ME2cDipY_max) :: dslist
  real(8), dimension (ME2cDipY_max) :: slist
  real(8), dimension (numorb_max,numorb_max) :: sm
  real(8), dimension (numorb_max,numorb_max) :: spm
  real(8), dimension (3,numorb_max,numorb_max) :: spmx
  logical switch
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
