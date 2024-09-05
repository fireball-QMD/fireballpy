subroutine cl_value (itype, cl)
  use iso_c_binding
  use M_system
  use M_fdata, only: nsshPP,lsshPP,num_orbPP,cl_PP
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
  if (index .ne. num_orbPP(itype)) then
    write (*,*) ' itype = ', itype
    write (*,*) ' index of orbitals for pseudopotential = ',index
    write (*,*) ' Program has num_orbPP = ', num_orbPP(itype)
    write (*,*) ' cl_value: index and num_orbPP DO NOT agree. '
    stop
  end if
  return
end
