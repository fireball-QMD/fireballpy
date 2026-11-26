subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_fdata, only: numXmax, numYmax, ME3c_max, nspecies
  implicit none
  integer, intent (in) :: index, iounit, isorp, maxtype, num_nonzero, numx, numy
  real(dp), intent (out), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3) :: xintegral
  integer :: ipoint, integral, jpoint

  do integral = 1, num_nonzero
    xintegral(numx+1:numXmax, :, integral, isorp, index) = 0.0_dp
    xintegral(1:numx, numy+1:numYmax, integral, isorp, index) = 0.0_dp
  end do

  do jpoint = 1, numy
    do ipoint = 1, numx
      read(iounit,*) (xintegral(ipoint,jpoint,integral,isorp,index), integral = 1, num_nonzero)
    end do
  end do
end subroutine readdata_3c

