subroutine assemble_F ()
  use iso_c_binding
  use M_system, only: natoms, neigh_j, neighn, nPP_j, nPPn, nPPx_j, nPPxn, fotnl, fanl, fotna, fana, faxc, faxc_ca, dxcdcc, ft, dusr, &
    & fotxc, fotxc_ca, faca, fotca, f3naa, f3nab, f3nac, f3nla, f3nlb, f3nlc, f3caa, f3cab, f3cac, flrew, f3xca_ca, f3xcb_ca, f3xcc_ca, &
    & f3xca, f3xcb, f3xcc, flrew_qmmm, fro, ftot, dxcv
  implicit none
  integer(c_long) iatom
  integer(c_long) ineigh
  integer(c_long) jatom
  integer(c_long) ix
  real(c_double) rms
  real(c_double) maximum
  real(c_double) minimum
  real(c_double), dimension (3, natoms) :: f3ca       ! three-center charged atom
  real(c_double), dimension (3, natoms) :: f3na       ! three-center neutral atom
  real(c_double), dimension (3, natoms) :: f3nl       ! three-center non-local
  real(c_double), dimension (3, natoms) :: f3xc       ! three-center xc
  real(c_double), dimension (3, natoms) :: f3xc_ca    ! three-center xc - charges
  real(c_double), dimension (3, natoms) :: fbs  ! total band-structure force
  real(c_double), dimension (3, natoms) :: fca  ! total charged atom force
  real(c_double), dimension (3, natoms) :: fcaatm     ! charged atom/atom force
  real(c_double), dimension (3, natoms) :: fcaot      ! charged atom/ontop force
  real(c_double), dimension (3, natoms) :: fna  ! total neutral atom force
  real(c_double), dimension (3, natoms) :: fnaatm     ! neutral atom/atom force
  real(c_double), dimension (3, natoms) :: fnaot      ! neutral atom/ontop force
  real(c_double), dimension (3, natoms) :: fnl  ! total non-local force
  real(c_double), dimension (3, natoms) :: fnlatm     ! non-local/atom force
  real(c_double), dimension (3, natoms) :: fnlot      ! non-local/ontop force
  real(c_double), dimension (3, natoms) :: fxc  ! total xc force
  real(c_double), dimension (3, natoms) :: fxc_ca     ! total xc force - charges
  real(c_double), dimension (3, natoms) :: fxcatm     ! xc atm/atm force
  real(c_double), dimension (3, natoms) :: fxcatm_ca  ! xc atm/atm force - charges
  real(c_double), dimension (3, natoms) :: fxcot      ! xc ontop force
  real(c_double), dimension (3, natoms) :: fxcot_ca   ! xc ontop force - charges
  do iatom = 1, natoms
    f3nl(:,iatom) = f3nla(:,iatom) + f3nlb(:,iatom) + f3nlc(:,iatom)
    f3na(:,iatom) = f3naa(:,iatom) + f3nab(:,iatom) + f3nac(:,iatom)
    f3ca(:,iatom) = f3caa(:,iatom) + f3cab(:,iatom) + f3cac(:,iatom)
    f3xc(:,iatom) = f3xca(:,iatom) + f3xcb(:,iatom) + f3xcc(:,iatom)
    f3xc_ca(:,iatom) = f3xca_ca(:,iatom) + f3xcb_ca(:,iatom) + f3xcc_ca(:,iatom)
  end do
  fcaatm = 0.0d0
  fnaatm = 0.0d0
  fcaot = 0.0d0
  fnaot = 0.0d0
  fxcatm = 0.0d0
  fxcatm_ca = 0.0d0
  fxcot = 0.0d0
  fxcot_ca = 0.0d0
  fnlatm = 0.0d0
  fnlot = 0.0d0
  dxcv = 0.0d0
  do iatom = 1, natoms
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      fnaatm(:,iatom) = fnaatm(:,iatom) + fana(:,ineigh,iatom)
      fnaatm(:,jatom) = fnaatm(:,jatom) - fana(:,ineigh,iatom)
      fcaatm(:,iatom) = fcaatm(:,iatom) + faca(:,ineigh,iatom)
      fcaatm(:,jatom) = fcaatm(:,jatom) - faca(:,ineigh,iatom)
      fxcatm(:,iatom) = fxcatm(:,iatom) + faxc(:,ineigh,iatom)
      fxcatm(:,jatom) = fxcatm(:,jatom) - faxc(:,ineigh,iatom)
      fxcatm_ca(:,iatom) = fxcatm_ca(:,iatom) + faxc_ca(:,ineigh,iatom)
      fxcatm_ca(:,jatom) = fxcatm_ca(:,jatom) - faxc_ca(:,ineigh,iatom)
      fnaot(:,iatom) = fnaot(:,iatom) + 2.0d0*fotna(:,ineigh,iatom)
      fnaot(:,jatom) = fnaot(:,jatom) - 2.0d0*fotna(:,ineigh,iatom)
      fcaot(:,iatom) = fcaot(:,iatom) + 2.0d0*fotca(:,ineigh,iatom)
      fcaot(:,jatom) = fcaot(:,jatom) - 2.0d0*fotca(:,ineigh,iatom)
      fxcot(:,iatom) = fxcot(:,iatom) + fotxc(:,ineigh,iatom)
      fxcot(:,jatom) = fxcot(:,jatom) - fotxc(:,ineigh,iatom)
      fxcot_ca(:,iatom) = fxcot_ca(:,iatom) + fotxc_ca(:,ineigh,iatom)
      fxcot_ca(:,jatom) = fxcot_ca(:,jatom) - fotxc_ca(:,ineigh,iatom)
      dxcv(:,iatom) = dxcv(:,iatom) + dxcdcc(:,ineigh,iatom)
      dxcv(:,jatom) = dxcv(:,jatom) - dxcdcc(:,ineigh,iatom)
    end do
  end do
  do iatom = 1, natoms
    do ineigh = 1, nPPn(iatom)
      jatom = nPP_j(ineigh,iatom)
      fnlatm(:,iatom) = fnlatm(:,iatom) + fanl(:,ineigh,iatom)
      fnlatm(:,jatom) = fnlatm(:,jatom) - fanl(:,ineigh,iatom)
     end do
  end do
  do iatom = 1, natoms
    do ineigh = 1, nPPxn(iatom)
      jatom = nPPx_j(ineigh,iatom)
      fnlot(:,iatom) = fnlot(:,iatom) + 2.0d0*fotnl(:,ineigh,iatom)
      fnlot(:,jatom) = fnlot(:,jatom) - 2.0d0*fotnl(:,ineigh,iatom)
    end do     ! end loop over PP-neighbors
  end do      ! end loop over atoms
  do iatom = 1, natoms
    fca(:,iatom) = f3ca(:,iatom) + fcaatm(:,iatom) + fcaot(:,iatom)
    fna(:,iatom) = f3na(:,iatom) + fnaatm(:,iatom) + fnaot(:,iatom)
    fnl(:,iatom) = f3nl(:,iatom) + fnlatm(:,iatom) + fnlot(:,iatom)
    fxc(:,iatom) = f3xc(:,iatom) + fxcatm(:,iatom) + fxcot(:,iatom)
    fxc_ca(:,iatom) = fxcatm_ca(:,iatom) + fxcot_ca(:,iatom)  + f3xc_ca(:,iatom)
  end do
  
  
  do iatom = 1, natoms
    fbs(:,iatom) = ft(:,iatom) + fna(:,iatom) + fnl(:,iatom) + fxc(:,iatom) + &
      & fca(:,iatom) + fxc_ca(:,iatom) + flrew(:,iatom) + flrew_qmmm(:,iatom)
  end do
  do iatom = 1, natoms
    ftot(:,iatom) = fbs(:,iatom) + dusr(:,iatom) + dxcv(:,iatom) + fro(:,iatom)
  end do
  rms = 0.0d0
  do iatom = 1, natoms
    do ix = 1,3
      rms = rms + ftot(ix,iatom)**2
    end do
  end do
  rms = sqrt(rms/(3*natoms))
  maximum = maxval(ftot)
  minimum = abs(minval(ftot))
  if (minimum .gt. maximum) maximum = minimum
 
end subroutine assemble_F
