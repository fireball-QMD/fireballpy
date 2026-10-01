subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmin, zmax, itype, in1, in2)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use math, only: math_interp_new
  use M_fdata, only: ME2c_max, nfofx, splineint_2c
  implicit none
  integer, intent (in) :: in1, in2, interaction, iounit, itype, num_nonzero, numz
  real(double), intent (in) :: zmin, zmax
  integer :: ipoint, integral
  real(double) :: dz
  real(double), allocatable :: xintegral_2c(:,:)
  allocate(xintegral_2c(num_nonzero, numz))
  do ipoint = 1, numz
    read (iounit,*) (xintegral_2c(integral,ipoint), integral = 1, num_nonzero)
  end do
  dz = (zmax - zmin)/real(numz - 1, double)
  do integral = 1, num_nonzero
    splineint_2c(integral,itype,in1,in2) = math_interp_new([(zmin + (ipoint - 1)*dz, ipoint = 1, numz)], &
      &                                                    xintegral_2c(integral,:), rval=0.0d0)
  end do
end subroutine readdata_2c
