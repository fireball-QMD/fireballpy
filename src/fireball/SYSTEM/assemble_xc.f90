subroutine assemble_xczw ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: idipole, iqmmm, natoms, eqmmm, Kscf, neigh_b, neigh_j, neighn, neigh_self, neighPP_self, neighPPn, &
    & neighPP_b, neighPP_j, ewaldqmmm, errno
  implicit none
  integer :: iatom, jatom, ineigh, mbeta, kforce
  if (Kscf .eq. 1) then

    call neighbors()
    call neighborsPP()
    call num_neigh_tot ()
    !call initneighbors ()
    call backnay ()
    call neighbors_pairs()
    call common_neighbors ()
    call common_neighborsPP ()
    if (errno .ne. 0) return
    kforce = 0
    call get_ewald (kforce) 
  end if ! end if (Kscf .eq. 1)

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
  !call assemble_zw_1c_na () !AQUI
  call assemble_xc_1c () !AQUI

  !  ----------------- assemble_2c ------------------------
  if (Kscf .eq. 1) then
    call assemble_sVNL ()
    call assemble_2c ()
    call assemble_2c_PP ()
  end if ! end if of Kscf = 1

  call average_rho() !AQUI

  !call average_ca_rho ()
  !call assemble_olsxc_on ()
  !call assemble_olsxc_off ()

  if (idipole .eq. 0) call assemble_ca_2c ()
  if (idipole .eq. 1) call assemble_ca_2c_dip ()
  

  call assemble_zw_on_na() !AQUI vxc = 0.0d0
  call assemble_zw_off_na() !AQUI
  call assemble_zw_2c_ct() !AQUI
  !-------------------- assemble_3c -------------------------
  if (Kscf .eq. 1) then
    call assemble_3c ()
    call assemble_3c_PP ()
    if (iqmmm .eq.1 ) then
      if (idipole .eq. 0) call assemble_qmmm ()
      if (idipole .eq. 1) call assemble_qmmm_dip ()
    else
      eqmmm = 0.0d0
      ewaldqmmm = 0.0d0
    end if
  end if
  if (idipole .eq. 0) call assemble_ca_3c ()
  if (idipole .eq. 1) call assemble_ca_3c_dip ()
  if (idipole .eq. 0) call assemble_lr ()
  if (idipole .eq. 1) call assemble_lr_dip ()


  call assemble_zw_3c_ct() !AQUI

  !Build H
  call buildh ()

end subroutine assemble_xczw
