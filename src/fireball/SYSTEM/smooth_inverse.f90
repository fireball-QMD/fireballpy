pure function sf(d)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: qmmm_rc1, qmmm_rc2, qmmm_width
  implicit none
  real(dp), intent(in) :: d
  real(dp) :: rc12, aux1, aux2, y, sf
  if (qmmm_rc1 < 0.0_dp) then
    aux1 = d*d
    aux1 = aux1*aux1
    aux2 = qmmm_rc1*qmmm_rc1
    aux2 = aux2*aux2
    sf = (aux1 - aux2)/(aux1*d + aux2*qmmm_rc1)
  else
    rc12 = qmmm_rc2 - qmmm_width
    if (d > qmmm_rc2) then
      sf = 0.0_dp
    else if (d > rc12) then
      y = (d - rc12) / qmmm_width
      sf = (1.0_dp - y*y*(3.0_dp - 2.0_dp*y))/d
    else if (d > qmmm_rc1) then
      sf = 1.0_dp/d
    else
      y = d/qmmm_rc1
      y = y*y
      sf = (3.28125_dp + y*(-5.46875_dp + y*(4.59375_dp + y*(-1.40625_dp))))/qmmm_rc1
    end if
  end if
end function
