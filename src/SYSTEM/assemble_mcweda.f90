subroutine assemble_mcweda ()
  use M_system
  implicit none
  integer iatom
  integer jatom
  integer ineigh
  integer mbeta
  integer kforce
  if (Kscf .eq. 1) then

    call neighbors()
    call neighborsPP()
    call num_neigh_tot ()
    call initneighbors ()
    call backnay ()
    call neighbors_pairs()
    call common_neighbors ()
    call common_neighborsPP ()

  end if ! end if (Kscf .eq. 1)
  kforce = 0
  call get_ewald (kforce) 

  neigh_self = -999
  do iatom = 1, natoms
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) neigh_self(iatom) = ineigh
    end do
  end do
  neighPP_self = -999
  do iatom = 1, natoms
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0)  neighPP_self(iatom) = ineigh
   end do
  end do

  ! assemble_1c
  call assemble_olsxc_1c ()

  !  ----------------- assemble_2c ------------------------
  if (Kscf .eq. 1) then
    call assemble_sVNL ()
    call assemble_2c ()
    call assemble_2c_PP ()
  end if ! end if of Kscf = 1

  call average_ca_rho ()
  call assemble_olsxc_on ()
  call assemble_olsxc_off ()

  if (idipole .eq. 0) call assemble_ca_2c ()
  if (idipole .eq. 1) call assemble_ca_2c_dip ()
  
  !-------------------- assemble_3c -------------------------
  if (Kscf .eq. 1) then
    call assemble_3c ()
    call assemble_3c_PP ()
    !if (iqmmm .eq.1 ) then
    !  if (idipole .eq. 0) call assemble_qmmm (nprocs, iordern)
    !  if (idipole .eq. 1) call assemble_qmmm_dip (nprocs, iordern)
    !else
    !  eqmmm = 0.0d0
    !  ewaldqmmm = 0.0d0
    !end if
  end if
  if (idipole .eq. 0) call assemble_ca_3c ()
  if (idipole .eq. 1) call assemble_ca_3c_dip ()
  if (idipole .eq. 0) call assemble_lr ()
  if (idipole .eq. 1) call assemble_lr_dip ()

  !Build H
  call buildh ()

end subroutine assemble_mcweda

