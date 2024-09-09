! The subroutine gets total number of neighbors (normal+PP) of each atom
subroutine num_neigh_tot ()
  use iso_c_binding
  use M_system, only: natoms, neigh_max, neigh_b, neigh_j, neighn, neighj_tot, neighb_tot, neighn_tot, numorb_max, neighPPn, neighPP_b, &
    & neighPP_j, neighj_aux, neighPP_max, num_neig_maxtot, hr_box
  implicit none
  integer(c_long),allocatable  :: neighb_aux(:,:)
  integer(c_long) :: iatom
  integer(c_long) :: jatom
  integer(c_long) :: jatomPP
  integer(c_long) :: num_neigh
  integer(c_long) :: num_neighPP
  integer(c_long) :: num_neig_tot
  integer(c_long) :: ineighPP
  integer(c_long) :: ineigh
  integer(c_long) :: count_neig
  integer(c_long) :: mbeta
  integer(c_long) :: mbetaPP

  if ( allocated (neighj_tot)) deallocate (neighj_tot)
  if ( allocated (neighb_tot)) deallocate (neighb_tot)
  if ( allocated (hr_box)) deallocate (hr_box)
  ! JOM-warning : these allocations seem arbitrary. we should improve
  allocate (neighj_aux(neigh_max+neighPP_max**2,natoms))
  allocate (neighb_aux(neigh_max+neighPP_max**2,natoms))
  neighj_aux = 0
  neighb_aux = 0

  do iatom = 1,natoms
    num_neigh = neighn(iatom)
    num_neighPP = neighPPn(iatom)
    num_neig_tot = num_neigh 
    neighj_aux(1:num_neigh,iatom) = neigh_j(1:num_neigh,iatom)
    neighb_aux(1:num_neigh,iatom) = neigh_b(1:num_neigh,iatom)
    do ineighPP = 1,num_neighPP
      count_neig = 0
      jatomPP = neighPP_j(ineighPP,iatom)
      mbetaPP = neighPP_b(ineighPP,iatom)
      do ineigh = 1, num_neigh
        jatom = neigh_j(ineigh,iatom)
        mbeta = neigh_b(ineigh,iatom)
        if ((jatomPP .eq. jatom .and. mbetaPP .eq. mbeta)) then
          count_neig = 1
        end if
      end do  ! ineigh 
      if (count_neig .eq. 0) then
        num_neig_tot = num_neig_tot + 1
        neighj_aux(num_neig_tot,iatom) = jatomPP 
        neighb_aux(num_neig_tot,iatom) = mbetaPP
      end if
    end do  ! ineighPP
    neighn_tot(iatom) = num_neig_tot
  end do  ! iatom
  num_neig_maxtot = maxval(neighn_tot(1:natoms))
  allocate (neighj_tot(num_neig_maxtot,natoms))
  allocate (neighb_tot(num_neig_maxtot,natoms))
  allocate (hr_box(numorb_max,numorb_max,natoms,0:num_neig_maxtot))
  neighj_tot = 0
  neighb_tot = 0
  hr_box = 0.0d0
  ! Loop over atoms
  do iatom = 1,natoms
    neighj_tot(1:neighn_tot(iatom),iatom) =  neighj_aux(1:neighn_tot(iatom),iatom)
    neighb_tot(1:neighn_tot(iatom),iatom) =  neighb_aux(1:neighn_tot(iatom),iatom)
  end do
  deallocate (neighb_aux)
  deallocate (neighj_aux)
  return
end subroutine num_neigh_tot
