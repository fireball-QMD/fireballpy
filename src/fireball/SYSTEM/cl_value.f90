subroutine cl_value (itype, cl)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max
  use M_fdata, only: nsshPP,lsshPP,cl_PP
  implicit none                  
  integer, intent (in) :: itype
  real(double), intent (out), dimension (numorb_max) :: cl
  integer imu
  integer index
  integer issh
  integer Lvalue
  integer Lmax
  cl(1:numorb_max) = 0.0d0
  index = 0
  do issh = 1, nsshPP(itype)
    Lvalue = lsshPP(issh,itype)
    Lmax = (2*Lvalue + 1)
    do imu = 1, Lmax
      index = index + 1
      cl(index) = cl_PP(issh,itype)
    end do
  end do
  return
end subroutine cl_value
