subroutine backnay ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, neigh_b, neigh_j, neighn, neigh_back, xl, errno
  implicit none
  integer :: iatom
  integer :: ineigh
  integer :: jatom
  integer :: jbeta
  integer :: jneigh
  integer :: katom
  integer :: mbeta
  integer :: iloop
  integer :: icount

  outer: do
    do iatom = 1, natoms
      inner: do ineigh = 1, neighn(iatom)
        mbeta = neigh_b(ineigh,iatom)
        jatom = neigh_j(ineigh,iatom)
        icount = 0
        do jneigh = 1, neighn(jatom)
          jbeta = neigh_b(jneigh,jatom)
          katom = neigh_j(jneigh,jatom)
          if (katom .eq. iatom .and. sum((xl(:,jbeta) + xl(:,mbeta))**2) .lt. 1.0d-5) then
            neigh_back (iatom,ineigh) = jneigh
            icount = icount + 1
          end if
        end do
        if (icount .eq. 1) cycle inner
        if (icount .gt. 1) then
          errno = 2
          return
        end if
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
