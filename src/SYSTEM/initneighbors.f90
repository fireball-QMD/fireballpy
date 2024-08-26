! This routine initializes the neighbors for the simulation. If the
! old neighbor file exists, then this will be used to determine the neighbors
! If not, then the neighbors will be calculated as normal. 
subroutine initneighbors ()
  use M_system
  implicit none
  integer :: iatom
  integer :: ineigh
  integer :: jatom
  integer :: jneigh
  integer :: jjneigh
  integer :: katom
  integer :: kneigh
  integer :: mbeta
  integer :: num_neigh
  integer :: num_neigh_vdw
  !   SELF nPP_self  
  ! What's the neighbor of the atom itself?
  nPP_self = -999
  do iatom = 1, natoms
    do ineigh = 1, nPPn(iatom)
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) nPP_self(iatom) = ineigh
    end do
  end do

  !   SELF nPPx_self  
  nPPx_self = -999
  do iatom = 1, natoms
    do ineigh = 1, nPPxn(iatom)
      mbeta = nPPx_b(ineigh,iatom)
      jatom = nPPx_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) nPPx_self(iatom) = ineigh
    end do
  end do

  return
end
