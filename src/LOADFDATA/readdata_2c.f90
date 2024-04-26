subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
  use M_fdata, only: ME2c_max, nfofx, xintegral_2c
  implicit none
  integer, intent (in) :: in1, in2, interaction, iounit, itype, num_nonzero, numz
  real*8, intent (in) :: zmax
  integer ipoint, integral
  real*8, dimension (ME2c_max, nfofx) :: gstore
  if (interaction .ne. 8) then
    do ipoint = 1, numz
      read (iounit,*) (gstore(integral,ipoint), integral = 1, num_nonzero)
    end do
    do ipoint = 1, numz
      do integral = 1, num_nonzero
        xintegral_2c(integral,ipoint,itype,in1,in2) = gstore(integral,ipoint)
      end do
    end do
  else
    do ipoint = 1, numz
      read (iounit,*) gstore(1,ipoint)
      xintegral_2c(1,ipoint,itype,in1,in2) = gstore(1,ipoint)
    end do
  end if
end subroutine readdata_2c
