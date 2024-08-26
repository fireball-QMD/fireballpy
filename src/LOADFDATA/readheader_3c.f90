subroutine readheader_3c (iounit, numx, numy, xmax, ymax)
  use M_constants, only: wp
  use M_fdata, only: numXmax, numYmax, errno3c
  implicit none
  integer, intent (in) :: iounit
  integer, intent (out) :: numx, numy
  real(wp), intent (out) :: xmax, ymax
  integer :: iline, nucZ1, nucZ2, nucZ3, nr, ntheta_in, nphi2
  real(wp) :: rc1a, rc2a, rc3a

  do iline = 1, 10
      read (iounit,*)
  end do
  read (iounit,*) nphi2, nr, ntheta_in
  read (iounit,*) ymax, numy
  read (iounit,*) xmax, numx
  read (iounit,*)
  read (iounit,*) nucZ1, rc1a
  read (iounit,*) nucZ2, rc2a
  read (iounit,*) nucZ3, rc3a
  read (iounit,*)
  if (numx .gt. numXmax .or. numy .gt. numYmax) then
    !write (*,*) ' Courseness too fine in 3c data files. '
    !write (*,*) ' numx = ', numx, ' numXmax = ', numXmax
    !write (*,*) ' numy = ', numy, ' numYmax = ', numYmax
    !write (*,*) ' Change numXmax and numYmax in MODULES/dimensions.f90! '
    errno3c = 1
    return
  end if
end subroutine readheader_3c
