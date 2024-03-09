integer function mpairnay (iatom, jatom, rdiff)
  use M_system
  implicit none
  integer, intent(in) :: iatom
  integer, intent(in) :: jatom
  real, intent(in), dimension (3) :: rdiff
  integer imatch
  integer ineigh
  integer jjatom
  integer mbeta
  real diff
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
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
  if (imatch .ne. 1) then
   write (*,*) ' imatch = ', imatch
   write (*,*) ' The variable imatch MUST be ONE! NO EXCEPTIONS '
   write (*,*) ' Bad imatch value in mpairnay.f90; must abort! '
   write (*,*) ' iatom = ', iatom, r1(:)
   write (*,*) ' jatom = ', jatom, rdiff(:) 
   stop
  end if
  if (mpairnay .le. 0) then
   write (*,*) ' Huh? The variable, mpairnay < 0, mpairnay = ', mpairnay
   write (*,*) ' It MUST be > 1. NO EXCEPTIONS '
   write (*,*) ' Bad mpairnay value in mpairnay.f90; must abort! '
   stop
  end if
  return
end function mpairnay

