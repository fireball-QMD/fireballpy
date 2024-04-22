subroutine doscentrosS (interaction, isub, iauxforce, in1, in2, in3, distance, eps, sx, spx)
  use M_system
  use M_fdata, only: index_maxS,nsh_max,MES_max,nssh
  implicit none
  integer, intent (in) :: iauxforce
  integer, intent (in) :: interaction
  integer, intent (in) :: isub
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: in3
  real(8), intent (inout) :: distance
  real(8), intent (in), dimension (3, 3) :: eps
  real(8), intent (out), dimension (nsh_max, nsh_max) :: sx
  real(8), intent (out), dimension (3, nsh_max, nsh_max) :: spx
  integer imu
  integer inu
  integer index
  real(8), dimension (3) :: eta
  real(8), dimension (MES_max) :: dslist
  real(8), dimension (MES_max) :: slist
  real(8), dimension (nsh_max,nsh_max) :: spm

  sx = 0.0d0
  if (iauxforce .eq. 1) spm = 0.0d0
  if (iauxforce .eq. 1) spx = 0.0d0
  do index = 1, index_maxS(in1,in3)
    if (interaction .ne. 20) then
      call interpolate_1d (interaction, isub, in1, in2, index, iauxforce, distance, slist(index), dslist(index))
    else
      call interpolate_1d (interaction, isub, in1, in3, index, iauxforce,  distance, slist(index), dslist(index))
    end if
  end do
  call recover_S (in1, in3, slist, sx)
  call recover_S (in1, in3, dslist, spm)
  if (iauxforce .eq. 1) then
    eta(:) = eps(:,3)
    do inu = 1, nssh(in3)
      do imu = 1, nssh(in1)
        if (distance .gt. 1.0d-3) then
          spx(:,imu,inu) = - eta(:)*spm(imu,inu)
        end if
      end do
    end do
  end if
  return
end subroutine doscentrosS
