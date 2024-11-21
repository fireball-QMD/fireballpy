subroutine backnay ()
  use iso_c_binding
  use M_system, only: natoms, neigh_b, neigh_j, neighn, neigh_back, xl
  implicit none
  integer(c_long) :: iatom
  integer(c_long) :: ineigh
  integer(c_long) :: jatom
  integer(c_long) :: jbeta
  integer(c_long) :: jneigh
  integer(c_long) :: katom
  integer(c_long) :: mbeta
  integer(c_long) :: iloop

  outer: do
    do iatom = 1, natoms
      inner: do ineigh = 1, neighn(iatom)
        mbeta = neigh_b(ineigh,iatom)
        jatom = neigh_j(ineigh,iatom)
        do jneigh = 1, neighn(jatom)
          jbeta = neigh_b(jneigh,jatom)
          katom = neigh_j(jneigh,jatom)
          if (katom .eq. iatom .and. sum((xl(:,jbeta) + xl(:,mbeta))**2) .lt. 1.0d-5) then
            neigh_back (iatom,ineigh) = jneigh
            cycle inner
          end if
        end do
        neighn(iatom) = neighn(iatom) - 1
        do iloop = ineigh, neighn(iatom)
          neigh_b(iloop,iatom) = neigh_b(iloop+1,iatom)
          neigh_j(iloop,iatom) = neigh_j(iloop+1,iatom)
        end do
        cycle outer
      end do inner
    end do
    exit outer
  end do outer
end subroutine backnay
