subroutine cl_value (itype, cl)
  use iso_c_binding
  use M_system, only: numorb_max
  use M_fdata, only: nsshPP,lsshPP,cl_PP
  implicit none                  
  integer(c_long), intent (in) :: itype
  real(c_double), intent (out), dimension (numorb_max) :: cl
  integer(c_long) imu
  integer(c_long) index
  integer(c_long) issh
  integer(c_long) Lvalue
  integer(c_long) Lmax
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
