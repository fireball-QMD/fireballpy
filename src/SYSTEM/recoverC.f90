! ===========================================================================
! This subroutine takes a 1D list of integrals and generates a 2x2 
! array with respect to the shells.
! ===========================================================================
subroutine recoverC (n1, n2, hlist, dhlist, hbox, dhbox)
  use M_fdata
  use M_system
  implicit none
  integer, intent (in) :: n1, n2
  real, intent(in), dimension (ME2c_max) :: hlist 
  real, intent(in), dimension (ME2c_max) :: dhlist 
  real, intent(out), dimension (nsh_max, nsh_max) :: hbox
  real, intent(out), dimension (nsh_max, nsh_max) :: dhbox

  integer index
  integer indexcoulomb
  integer ii
  integer kk
 
  ! The algorithm is simply to go along the vector, and fill in the first row 
  ! of a matrix. Once the nssh(2) elements are absorbed, we increment the row 
  ! index ii by one, and reset the column index kk back to one.
  ii = 1
  kk = 0

  ! Loop over all the non-zero integrals for this interaction:
  ! n1 = number of shells on atom 1, and n2 = number of shells on atom2.
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
 
