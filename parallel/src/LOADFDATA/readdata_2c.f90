subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
  use M_fdata, only: ME2c_max, nfofx
  implicit none
  integer, intent (in) :: in1, in2, interaction, iounit, itype, num_nonzero, numz
  real*8, intent (in) :: zmax
  integer :: ipoint, integral
  real*8, dimension (ME2c_max, nfofx) :: xintegral_2c
  !if (interaction .ne. 8) then
  xintegral_2c = 0.0d0
  do ipoint = 1, numz
    read (iounit,*) (xintegral_2c(integral,ipoint), integral = 1, num_nonzero)
  end do
  do integral = 1, num_nonzero
    call buildspline_1d(integral, numz, itype, in1, in2, zmax, interaction, xintegral_2c)
  end do
  !else
  !  do ipoint = 1, numz
  !    read (iounit,*) xintegral_2c(1,ipoint,itype,in1,in2)
  !  end do
  !end if
end subroutine readdata_2c
