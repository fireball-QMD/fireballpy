subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
  use M_fdata, only: ME2c_max, xintegral_2c
  implicit none
  integer, intent (in) :: in1, in2, interaction, iounit, itype, num_nonzero, numz
  real*8, intent (in) :: zmax
  integer ipoint, integral
  real*8, dimension (ME2c_max) :: gstore
  if (interaction .ne. 8) then
    do ipoint = 1, numz
      read (iounit,*) (gstore(integral), integral = 1, num_nonzero)
      do integral = 1, num_nonzero
        xintegral_2c(integral,ipoint,itype,in1,in2) = gstore(integral)
      end do
    end do
  else
    do ipoint = 1, numz
      read (iounit,*) xintegral_2c(1,ipoint,itype,in1,in2)
    end do
  end if
end subroutine readdata_2c
