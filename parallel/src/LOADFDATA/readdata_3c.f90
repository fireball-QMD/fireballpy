subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
  use M_fdata, only: numXmax, numYmax, ME3c_max, nspecies
  implicit none
  integer, intent (in) :: index, iounit, isorp, maxtype, num_nonzero, numx, numy
  real*8, intent (out), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies*nspecies*nspecies) :: xintegral
  integer ipoint, integral, jpoint
  real*8, dimension (ME3c_max, numXmax, numYmax) :: gstore

  do jpoint = 1, numy
    do ipoint = 1, numx
    read (iounit,*) (gstore(integral,ipoint,jpoint), integral = 1, num_nonzero)
    end do
  end do
  do jpoint = 1, numy
    do ipoint = 1, numx
    do integral = 1, num_nonzero
      xintegral(ipoint,jpoint,integral,isorp,index) = gstore(integral,ipoint,jpoint)
    end do
    end do
  end do
end subroutine readdata_3c

