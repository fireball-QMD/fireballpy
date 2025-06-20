! This subroutine takes a 1D list of integrals and generates a 2x2 array with respect to the shells.
subroutine recoverC (n1, n2, hlist, dhlist, hbox, dhbox)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: ME2c_max, nsh_max
  implicit none
  integer, intent (in) :: n1, n2
  real(double), intent(in), dimension (ME2c_max) :: hlist 
  real(double), intent(in), dimension (ME2c_max) :: dhlist 
  real(double), intent(out), dimension (nsh_max, nsh_max) :: hbox
  real(double), intent(out), dimension (nsh_max, nsh_max) :: dhbox

  integer index
  integer indexcoulomb
  integer ii
  integer kk
 
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
end subroutine recoverC
