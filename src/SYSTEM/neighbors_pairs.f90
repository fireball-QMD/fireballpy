subroutine neighbors_pairs ()
  use M_system
  implicit none
  integer :: num_pairs
  integer :: iatom
  integer :: ineigh
  integer :: jatom
  integer :: jneigh
  num_pairs = 0
  if (icluster .eq. 1) then
    do iatom = 1,natoms
      do ineigh = 1,neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        if (jatom .gt. iatom) then   !Doubt.. .gt. or .ge.? what's best?
          num_pairs = num_pairs+1
          jneigh = neigh_back(iatom,ineigh)
          neigh_pair_a1(num_pairs) = iatom
          neigh_pair_a2(num_pairs) = jatom
          neigh_pair_n1(num_pairs) = ineigh
          neigh_pair_n2(num_pairs) = jneigh
        end if !if jatom .gt. iatom
      end do !end do ineigh = 1,neighn(iatom) 
    end do !end do iatom = 1,natoms
  else  !if icluster .eq. 1
  end if !end if icluster .eq. 1
  tot_pairs = num_pairs    !tot_pairs stores the total number of non-repeated
  return 
end
