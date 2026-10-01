subroutine interpolate_1d (interaction, isub, in1, in2, non2c, ioption, xin, yout, dfdx)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: ind2c, numz2c, z2cmax, splineint_2c
  implicit none
  integer, intent(in) :: interaction, isub, in1, in2, non2c, ioption ! Derivative or not
  real(double), intent(in)  :: xin
  real(double), intent(out) :: yout, dfdx
  integer :: jxx
  jxx = ind2c(interaction,isub)
  yout = splineint_2c(non2c,jxx,in1,in2)%f(xin)
  if (ioption == 1) dfdx = splineint_2c(non2c,jxx,in1,in2)%f(xin, order=1)
end subroutine interpolate_1d
