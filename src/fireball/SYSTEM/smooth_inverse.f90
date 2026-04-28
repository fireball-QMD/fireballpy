pure function sf(d)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: qmmm_rc1, qmmm_rc2, qmmm_width
  implicit none
  real(double), intent(in) :: d
  real(double) :: rc12, aux1, aux2, y, sf
  if (qmmm_rc1 < 0.0d0) then
    aux1 = d*d
    aux1 = aux1*aux1
    aux2 = qmmm_rc1*qmmm_rc1
    aux2 = aux2*aux2
    sf = (aux1 - aux2)/(aux1*d + aux2*qmmm_rc1)
  else
    rc12 = qmmm_rc2 - qmmm_width
    if (d > qmmm_rc2) then
      sf = 0.0d0
    else if (d > rc12) then
      y = (d - rc12) / qmmm_width
      sf = (1.0d0 - y*y*(3.0d0 - 2.0d0*y))/d
    else if (d > qmmm_rc1) then
      sf = 1.0d0/d
    else
      y = d/qmmm_rc1
      y = y*y
      sf = (3.28125d0 + y*(-5.46875d0 + y*(4.59375d0 + y*(-1.40625d0))))/qmmm_rc1
    end if
  end if
end function
