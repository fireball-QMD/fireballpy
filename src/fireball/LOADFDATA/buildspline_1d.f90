subroutine buildspline_1d (integral, numz, itype, in1, in2, zmax, interaction, xintegral_2c)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: ind2c, ME2c_max, nfofx, splineint_2c
  implicit none
  integer, intent(in) :: integral, numz, itype, in1, in2, interaction
  real(double), intent(in) :: zmax
  real(double), dimension(ME2c_max,nfofx), intent(in) :: xintegral_2c
  integer :: iorder, numz_used
  real(double) :: h, hi, zmin, w
  real(double), dimension (numz) :: b, d, z

  zmin = 0.0d0
  h = (zmax - zmin)/(numz - 1)
  hi = 1.0d0/h
  numz_used = numz
  do iorder = numz-1,3,-1
    if (xintegral_2c(integral,iorder) .ne. 0) then
      numz_used = iorder + 1
      exit
    end if
  end do

  b(1) = 2.0d0
  d(1) = 0.0d0
  do iorder=2,numz_used-1
    w = 1.0d0/b(iorder-1)
    b(iorder) = 4.0d0 - w
    d(iorder) = 3.0d0*(xintegral_2c(integral,iorder+1) - 2.0d0*xintegral_2c(integral,iorder) + xintegral_2c(integral,iorder-1)) - w*d(iorder-1)
  end do
  w = 1.0d0/b(numz_used-1)
  b(numz_used) = 2.0d0 - w
  ! Slight correction for 1/r derivative JOM-ENRIQUE
  if(((interaction .eq. 12) .or. (interaction .eq. 4)) .and. (itype .ne. ind2c(4,0))) then
    d(numz_used) = -3.0d0*h*(xintegral_2c(integral,numz_used) - xintegral_2c(integral,numz_used-1))/zmax - w*d(numz_used-1)
  else
    d(numz_used) = -w*d(numz_used-1)
  end if

  z(numz_used) = d(numz_used)/b(numz_used)
  do iorder=numz_used-1,1,-1
    z(iorder) = (d(iorder) - z(iorder+1))/b(iorder)
  end do

  do iorder=1,numz_used
    splineint_2c(1,integral,iorder,itype,in1,in2) = xintegral_2c(integral,iorder)
    splineint_2c(2,integral,iorder,itype,in1,in2) = xintegral_2c(integral,iorder+1) - xintegral_2c(integral,iorder) &
      & - 0.333333333333333d0*(z(iorder+1) + 2.0*z(iorder))
    splineint_2c(3,integral,iorder,itype,in1,in2) = z(iorder)
    splineint_2c(4,integral,iorder,itype,in1,in2) = 0.333333333333333d0*(z(iorder+1) - z(iorder))
  end do
end subroutine buildspline_1d
