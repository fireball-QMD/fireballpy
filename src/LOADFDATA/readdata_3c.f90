subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
  use M_fdata, only: numXmax, numYmax, ME3c_max, nspecies
  implicit none
  integer, intent (in) :: index, iounit, isorp, maxtype, num_nonzero, numx, numy
  real*8, intent (out), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3) :: xintegral
  integer :: ipoint, integral, jpoint

  do jpoint = 1, numy
    do ipoint = 1, numx
      read(iounit,*) (xintegral(ipoint,jpoint,integral,isorp,index), integral = 1, num_nonzero)
    end do
  end do
end subroutine readdata_3c

