subroutine allocate_system ()
  use M_system 
  use M_constants
  use M_fdata, only: nssh, rcutoff, rc_PP, nspecies, symbolA ,isorpmax
  use M_fdata, only: num_orb, Qneutral, lssh, nsshPP, lsshPP,  nsh_max
  use M_fdata, only: ME3c_max,ideriv_max,numXmax,numYmax
  use M_fdata, only: numy3c_xc3c
  implicit none
  integer:: iatom
  integer:: jatom
  integer:: matom
  integer:: ineigh
  integer:: i,iorb
  integer:: mbeta
  integer:: num_neigh
  integer:: in1
  integer:: imu, qmu
  integer:: issh
  integer:: in2
  integer:: ispec
  integer:: numorb
  integer:: numorbPP_max
  real :: rcutoff_i
  real :: rcutoff_j
  real:: rcutoff_
  real:: distance2
  real:: range2

  if (.not. allocated (ratom)) allocate (ratom (3, natoms))
  if (.not. allocated (imass)) allocate (imass (natoms))
  if (.not. allocated (symbol)) allocate (symbol (natoms))

  allocate (degelec (natoms))
  allocate (nelectron(natoms))
  allocate (Qin(nsh_max, natoms))
  allocate (Qinmixer(nsh_max*natoms))
  allocate (Qoutmixer(nsh_max*natoms))
  allocate (QLowdin_TOT (natoms))
  allocate (QMulliken_TOT (natoms))
  allocate (Qout(nsh_max, natoms))
  allocate (dq(nspecies))
  allocate (Q0_TOT(natoms))
  allocate (Q_partial(natoms))
  allocate (dq_DP(natoms))

  call initboxes ()

  ! Count the orbitals
  norbitals = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    norbitals = norbitals + num_orb(in1)
  end do

  ! Count total n of shells in the system.
  nssh_tot = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      nssh_tot = nssh_tot + 1
    end do
  end do

  numorb_max = 0
  do in1 = 1, nspecies
    numorb = 0
    do issh = 1, nssh(in1)
      numorb = numorb + 2*lssh(issh,in1) + 1
    end do
    if (numorb .gt. numorb_max) numorb_max = numorb
  end do


  ! Qneutral_total
  do iatom = 1, natoms
    Q0_TOT(iatom) = 0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
    end do
  end do


  ! By default the input charges are initialized to the neutral atom charges
  Qin = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Qin(issh,iatom) = Qneutral(issh,in1)
    end do
  end do


  ztot = 0.0d0
  nelectron = 0.0d0
  do iatom = 1, natoms
   in1 = imass(iatom)
   do issh = 1, nssh(in1)
    ztot = ztot + Qneutral(issh,in1)
    nelectron(iatom) = nelectron(iatom) + Qneutral(issh,in1)
   end do
  end do

  ! Calculate degelec. 
  degelec(1) = 0
  do iatom = 2, natoms
    degelec(iatom) = 0
    in1 = imass(iatom - 1)
    degelec(iatom) = degelec(iatom - 1) + num_orb(in1)
  end do
  numorbPP_max = 0
  do in1 = 1, nspecies
    numorb = 0
    do issh = 1, nsshPP(in1)
      numorb = numorb + 2*lsshPP(issh,in1) + 1
    end do
  if (numorb .gt.  numorbPP_max) numorbPP_max = numorb
  end do
  if (numorbPP_max .gt.  numorb_max) numorb_max = numorbPP_max

  allocate (getmssh(norbitals))
  allocate (getlssh(norbitals))
  allocate (getissh(norbitals))
  allocate (getiatom(norbitals))
 

             
  imu=0
  do iatom=1,natoms
    do issh=1,nssh(imass(iatom))
      do iorb = 1, 2*lssh(issh,imass(iatom))+1 
        imu=imu+1
        getissh(imu)=issh
        getlssh(imu)=lssh(issh,imass(iatom))
        getiatom(imu)=iatom
        getmssh(imu)=iorb 
      end do
   end do
  end do

 
  if (icluster .eq. 1) mbeta_max = 0
  neigh_max = -99
  do iatom = 1, natoms
    num_neigh = 0
    rcutoff_i = 0.0d0
    in1 = imass(iatom)
    do imu = 1, nssh(in1)
      if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do
    do mbeta = 0, mbeta_max
      do jatom = 1, natoms
        rcutoff_j = 0.0d0
        in2 = imass(jatom)
        do imu = 1, nssh(in2)
          if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
        end do
        distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2 + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2  + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
        if (distance2 .le. range2) num_neigh = num_neigh + 1
      end do
    end do
  end do

  neigh_max = max(neigh_max, num_neigh)

 
  ! find_neighPP_max
  if (icluster .eq. 1) mbeta_max = 0
  neighPP_max = -99
  do iatom = 1, natoms
    num_neigh = 0
    rcutoff_i = 0.0d0
    in1 = imass(iatom)
    do imu = 1, nssh(in1)
      if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do

    ! Loop over all possible neighbors (VNL atoms)
    do mbeta = 0, mbeta_max
      do jatom = 1, natoms
        in2 = imass(jatom)
        rcutoff_j = rc_PP(in2)
        distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
        &         + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2  &
        &         + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2

        ! Add a small displacement to the sum of cutoffs. 
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2

        if (distance2 .le. range2) then
          num_neigh = num_neigh + 1
        end if
      end do
    end do
    neighPP_max = max(neighPP_max, num_neigh)
  end do

  do iatom = 1, natoms
    num_neigh = 0
    in1 = imass(iatom)
    rcutoff_i = rc_PP(in1)
    ! Loop over all possible neighbors (VNL atoms)
    do mbeta = 0, mbeta_max
      do jatom = 1, natoms
        rcutoff_j = 0.0d0
        in2 = imass(jatom)
        do imu = 1, nssh(in2)
          if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
        end do
        distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
        &         + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2  &
        &         + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
        if (distance2 .le. range2) then
          num_neigh = num_neigh + 1
        end if
      end do
    end do
    neighPP_max = max(neighPP_max, num_neigh)
  end do

  allocate (neigh_b (neigh_max, natoms))
  allocate (neigh_j (neigh_max, natoms))
  allocate (neigh_comb (2, neigh_max**2, natoms))
  allocate (neigh_comj (2, neigh_max**2, natoms))
  allocate (neigh_comm (neigh_max**2, natoms))
  allocate (neigh_comn (natoms))
  allocate (neigh_back (natoms, neigh_max))
  allocate (neigh_self (natoms))
  allocate (nPP_b (neighPP_max, natoms))
  allocate (nPP_j (neighPP_max, natoms))
  allocate (nPP_map (neighPP_max, natoms))
  allocate (nPPn (natoms))
  allocate (nPP_self (natoms))
  allocate (nPPx_b (neighPP_max, natoms))
  allocate (nPPx_j (neighPP_max, natoms))
  allocate (nPPx_map (neighPP_max, natoms))
  allocate (nPPx_point (neighPP_max, natoms))
  allocate (nPPxn (natoms))
  allocate (neigh_pair_a1 (neigh_max*natoms))
  allocate (neigh_pair_a2 (neigh_max*natoms))
  allocate (neigh_pair_n1 (neigh_max*natoms))
  allocate (neigh_pair_n2 (neigh_max*natoms))
  allocate (nPPx_self (natoms))
  allocate (vxc_1c (numorb_max, numorb_max, neigh_max, natoms))
  allocate (neighn (natoms))
  allocate (neighn_tot (natoms))

! neighPP
  allocate (neighPP_b (neighPP_max**2, natoms))
  allocate (neighPP_j (neighPP_max**2, natoms))
  allocate (neighPPn (natoms))
! 3. party PPcommon pairs
  allocate (neighPP_comb (2, neighPP_max**2, natoms))
  allocate (neighPP_comj (2, neighPP_max**2, natoms))
  allocate (neighPP_comm (neighPP_max**2, natoms))
  allocate (neighPP_comn (natoms))
  allocate (neighPP_self (natoms))

! Total neighbor list (mapping together neigh and neighPP part)
  allocate (neighj_tot (neigh_max+neighPP_max, natoms))
  allocate (neighb_tot (  neigh_max+neighPP_max, natoms))

  
 
  allocate(ioccupy(norbitals))
  allocate(ioccupy_k (norbitals, nkpoints))
  allocate(foccupy (norbitals, nkpoints))
  allocate (eigen_k (norbitals, nkpoints))

  
  allocate (dewald (3, natoms, natoms))
  allocate (fewald (3, natoms))
  allocate (ewald (natoms, natoms))


  allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
  allocate (t_mat (numorb_max, numorb_max, neigh_max, natoms))
  allocate (h_mat (numorb_max, numorb_max, neigh_max, natoms))
  allocate (sp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
  allocate (tp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
  allocate (dipcm (3, numorb_max, numorb_max))
  allocate (dipc (3, numorb_max, numorb_max, neigh_max, natoms))
  allocate (dip (numorb_max, numorb_max, neigh_max, natoms)) 
  allocate (dippcm (3, 3, numorb_max, numorb_max))
  allocate (dippc (3, 3, numorb_max, numorb_max, neigh_max, natoms))
  allocate (sm_mat (nsh_max, nsh_max, neigh_max, natoms))
  allocate (ewaldqmmm (numorb_max, numorb_max, neigh_max,natoms))
  allocate (cape (numorb_max, numorb_max, neigh_max, natoms))

  allocate (rho (numorb_max, numorb_max, neigh_max, natoms))
  allocate (rhoPP (numorb_max, numorb_max, neighPP_max**2, natoms))
  allocate (dusr (3, natoms))
  

  allocate (bbnkre_o(norbitals,norbitals,nkpoints))
  allocate (arho_off (nsh_max, nsh_max, neigh_max, natoms))
  allocate (arhoij_off (nsh_max, nsh_max, neigh_max, natoms))
  allocate (rho_off (numorb_max, numorb_max, neigh_max, natoms))
  allocate (rhoij_off (numorb_max, numorb_max, neigh_max, natoms))
  allocate (arho_on (nsh_max, nsh_max, natoms))
  allocate (arhoi_on (nsh_max, nsh_max, natoms))
  allocate (rho_on (numorb_max, numorb_max, natoms))
  allocate (rhoi_on (numorb_max, numorb_max, natoms))
  allocate (vxc_ca (numorb_max, numorb_max, neigh_max, natoms))
  allocate (vxc (numorb_max, numorb_max, neigh_max, natoms))
  allocate (vca (numorb_max, numorb_max, neigh_max, natoms))
  allocate (ewaldsr (numorb_max, numorb_max, neigh_max, natoms))
  allocate (ewaldlr (numorb_max, numorb_max, neigh_max, natoms))

  allocate (spm_mat (3, nsh_max, nsh_max, neigh_max, natoms))
  allocate (arhop_off (3, nsh_max, nsh_max, neigh_max, natoms))
  allocate (arhopij_off (3, nsh_max, nsh_max, neigh_max, natoms)) 
  allocate (rhop_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
  allocate (rhopij_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
  allocate (arhop_on (3, nsh_max, nsh_max, neigh_max, natoms))
  allocate (rhop_on (3, numorb_max, numorb_max, neigh_max, natoms)) 
  ! OLSXC double count corr forces
  allocate (dxcdcc (3, neigh_max, natoms))

  allocate (vna (numorb_max, numorb_max, neigh_max, natoms))
  allocate (vnl (numorb_max, numorb_max, neighPP_max**2, natoms))
  allocate (sVNL (numorb_max, numorb_max, neighPP_max, natoms))
  allocate (spVNL (3, numorb_max, numorb_max, neighPP_max, natoms))

  call initamat()
  
  

  scf_achieved = .false.
  delk = 0.0d0
  delk(1,1) = 1.0d0
  delk(2,2) = 1.0d0
  delk(3,3) = 1.0d0
  xlevi = 0.0d0
  xlevi(1,2,3) = 1.0d0
  xlevi(1,3,2) = -1.0d0
  xlevi(3,1,2) = 1.0d0
  xlevi(3,2,1) = -1.0d0
  xlevi(2,3,1) = 1.0d0
  xlevi(2,1,3) = -1.0d0

  allocate (fana (3, neigh_max, natoms))
  allocate (fotna (3, neigh_max, natoms)) 
  allocate (ft (3, natoms))
  allocate (fanl (3, neighPP_max, natoms))
  allocate (fotnl (3, neighPP_max, natoms))
  allocate (faxc_ca (3, neigh_max, natoms)) 
  allocate (faxc (3, neigh_max, natoms))
  allocate (fotxc_ca (3, neigh_max, natoms))
  allocate (fotxc (3, neigh_max, natoms))
  allocate (faca (3, neigh_max, natoms))  
  allocate (fotca (3, neigh_max, natoms))
  allocate (f3naa (3, natoms))  
  allocate (f3nab (3, natoms))  
  allocate (f3nac (3, natoms))
  allocate (f3nla (3, natoms))  
  allocate (f3nlb (3, natoms))  
  allocate (f3nlc (3, natoms))
  allocate (f3caa (3, natoms))
  allocate (f3cab (3, natoms))
  allocate (f3cac (3, natoms))
  allocate (flrew (3, natoms))
  allocate (f3xca (3, natoms))
  allocate (f3xcc (3, natoms))
  allocate (f3xcb (3, natoms))
  allocate (f3xca_ca (3, natoms))
  allocate (f3xcb_ca (3, natoms))
  allocate (f3xcc_ca (3, natoms))
  allocate (flrew_qmmm (3, natoms))
  allocate (fro (3, natoms))
  allocate (ftot (3, natoms))
  allocate (dxcv (3, natoms))  
end subroutine
