subroutine readheader_3c (iounit, numx, numy, xmax, ymax)
  use iso_c_binding
  use M_fdata, only: numXmax, numYmax
  implicit none
  integer(c_long), intent (in) :: iounit
  integer(c_long), intent (out) :: numx, numy
  real(c_double), intent (out) :: xmax, ymax
  integer(c_long) :: iline, nucZ1, nucZ2, nucZ3, nr, ntheta_in, nphi2
  real(c_double) :: rc1a, rc2a, rc3a

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
end subroutine readheader_3c
