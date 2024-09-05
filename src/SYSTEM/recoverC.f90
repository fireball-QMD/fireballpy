! This subroutine takes a 1D list of integrals and generates a 2x2 array with respect to the shells.
subroutine recoverC (n1, n2, hlist, dhlist, hbox, dhbox)
  use iso_c_binding
  use M_fdata
  use M_system
  implicit none
  integer(c_long), intent (in) :: n1, n2
  real(c_double), intent(in), dimension (ME2c_max) :: hlist 
  real(c_double), intent(in), dimension (ME2c_max) :: dhlist 
  real(c_double), intent(out), dimension (nsh_max, nsh_max) :: hbox
  real(c_double), intent(out), dimension (nsh_max, nsh_max) :: dhbox

  integer(c_long) index
  integer(c_long) indexcoulomb
  integer(c_long) ii
  integer(c_long) kk
 
  ii = 1
  kk = 0

  indexcoulomb = n1*n2
  do index = 1, indexcoulomb
    kk = kk + 1
    hbox(ii,kk) = hlist(index)
    dhbox(ii,kk) = dhlist(index)
    if (mod(index,n2) .eq. 0) then
      ii = ii + 1
      kk = kk - n2
    end if
  end do
  return
end subroutine
 
