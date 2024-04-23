subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
  use M_fdata
  implicit none
  integer, intent (in) :: index
  integer, intent (in) :: iounit
  integer, intent (in) :: isorp
  integer, intent (in) :: maxtype
  integer, intent (in) :: num_nonzero
  integer, intent (in) :: numx, numy 
  real*8, intent (inout), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3) :: xintegral
  integer ipoint
  integer integral
  integer jpoint
  real*8, dimension (ME3c_max, numXmax, numYmax) :: gstore
  real*8, allocatable, save, dimension(:,:) :: binomial
  integer maxmax
  real*8, external :: factorial
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
  return
end subroutine readdata_3c

