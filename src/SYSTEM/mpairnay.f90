integer(c_long) function mpairnay (iatom, jatom, rdiff)
  use iso_c_binding
  use M_system, only: ratom, nPP_b, nPP_j, nPPn, xl
  implicit none
  integer(c_long), intent(in) :: iatom
  integer(c_long), intent(in) :: jatom
  real(c_double), intent(in), dimension (3) :: rdiff
  integer(c_long) imatch
  integer(c_long) ineigh
  integer(c_long) jjatom
  integer(c_long) mbeta
  real(c_double) diff
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r21
  mpairnay = 0
  r1(:) = ratom(:,iatom)
  imatch = 0
  do ineigh = 1, nPPn(iatom)     ! <==== loop 2 over i's neighbors
   jjatom = nPP_j(ineigh,iatom)
   if (jjatom .eq. jatom) then
    mbeta = nPP_b(ineigh,iatom)
    r2(:) = ratom(:,jatom) + xl(:,mbeta)
    r21(:) = r2(:) - r1(:)
    diff = (r21(1) - rdiff(1))**2 + (r21(2) - rdiff(2))**2       &
    + (r21(3) - rdiff(3))**2
    if (diff .lt. 0.0001d0) then
     imatch = imatch + 1
     mpairnay = ineigh
    end if
   end if
  end do
  return
end function mpairnay
