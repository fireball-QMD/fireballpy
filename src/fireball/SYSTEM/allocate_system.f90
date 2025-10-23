subroutine allocate_system ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: xlevi, delk
  use M_system, only: icluster, natoms, ratom, degelec, imass, mbeta_max, neigh_max, getmssh, getlssh, getissh, vxc_1c, &
    & Q0_TOT, nelectron, eigen_k, norbitals, nkpoints, getiatom, ioccupy_k, ioccupy, foccupy, cape, rhoPP, ztot, nssh_tot, ewald, &
    & dewald, fewald, ewaldsr, dip, dipp, neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, neigh_comm, neigh_comn, neigh_back, &
    & neigh_self, nPP_b, nPP_j, nPP_map, nPPn, nPP_self, nPPx_b, nPPx_j, nPPx_map, nPPx_point, nPPxn, nPPx_self, neigh_pair_a1, &
    & neigh_pair_a2, neigh_pair_n1, neigh_pair_n2, neighj_tot, neighb_tot, neighn_tot, numorb_max, neighPP_self, neighPPn, neighPP_b, &
    & neighPP_j, sVNL, spVNL, sp_mat, tp_mat, dipcm, dippcm, dippc, vnl, neighPP_comn, neighPP_comm, neighPP_comj, neighPP_comb, &
    & neighPP_max, Qin, Qinmixer, Qout, Qoutmixer, dq, Q_partial, QLowdin_TOT, QMulliken_TOT, dq_DP, vxc, vxc_ca, rho, rho_off, &
    & rhoij_off, s_mat, sm_mat, spm_mat, rho_on, arho_on, rhoi_on, arhoi_on, arhop_on, rhop_on, arhoij_off, arho_off, arhopij_off, &
    & arhop_off, rhop_off, rhopij_off, vca, ewaldlr, h_mat, t_mat, vna, ewaldqmmm, dipc, xl, fotnl, fanl, fotna, fana, faxc, faxc_ca, &
    & dxcdcc, ft, dusr, fotxc, fotxc_ca, faca, fotca, f3naa, f3nab, f3nac, f3nla, f3nlb, f3nlc, f3caa, f3cab, f3cac, flrew, f3xca_ca, &
    & f3xcb_ca, f3xcc_ca, f3xca, f3xcb, f3xcc, flrew_qmmm, fro, ftot, dxcv, norbitals_new, qstate, bbnkre, bbnkim, igamma, &
    & g_h, g_xc, f_xc, exc_aa, vxc_aa, get_orb_ofshell, get_l_ofshell, get_issh_ofshell, get_iatom_ofshell, get_shell_oforb, &
    & g_h_shell, g_xc_shell, f_xc_shell, exc_aa_shell, vxc_aa_shell
  use M_fdata, only: nssh, rcutoff, rc_PP, nspecies, num_orb, Qneutral, lssh, nsshPP, lsshPP,  nsh_max, numXmax, numYmax
!  use M_fdata, only: numy3c_xc3c, ideriv_max
  implicit none
  integer:: iatom
  integer:: jatom
  integer:: iorb
  integer:: mbeta
  integer:: num_neigh
  integer:: in1
  integer:: imu
  integer:: alpha
  integer:: issh
  integer:: in2
  integer:: numorb
  integer:: numorbPP_max
  real(double) :: rcutoff_i
  real(double) :: rcutoff_j
  real(double) :: distance2
  real(double) :: range2

  !Shift the coordinates none of the atoms fall on (0.0, 0.0, 0.0)
  !shifter(1) = 4.0d0*atan(1.0d0)    ! pi
  !shifter(2) = 1.0/exp(1.0d0)       ! 1/e
  !shifter(3) = sqrt(2.0d0)          ! square root of 2
  !ishiftO = 0
  !do iatom = 1, natoms
  !  distance = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2
  !  if (distance .lt. 1.0d-8) then
  !    ishiftO = 1
  !    exit
  !  end if
  !end do

  !if (ishiftO .eq. 1) then
  !  do iatom = 1, natoms
  !    ratom(:,iatom) = ratom(:,iatom) + shifter
  !  end do
  !end if
 

  if (allocated(degelec)) deallocate(degelec)
  allocate (degelec (natoms))
  if (allocated(nelectron)) deallocate(nelectron)
  allocate (nelectron(natoms))
  if (allocated(Qin)) deallocate(Qin)
  allocate (Qin(nsh_max, natoms))
  if (allocated(Qinmixer)) deallocate(Qinmixer)
  allocate (Qinmixer(nsh_max*natoms))
  if (allocated(Qoutmixer)) deallocate(Qoutmixer)
  allocate (Qoutmixer(nsh_max*natoms))
  if (allocated(QLowdin_TOT)) deallocate(QLowdin_TOT)
  allocate (QLowdin_TOT (natoms))
  if (allocated(QMulliken_TOT)) deallocate(QMulliken_TOT)
  allocate (QMulliken_TOT (natoms))
  if (allocated(Qout)) deallocate(Qout)
  allocate (Qout(nsh_max, natoms))
  if (allocated(dq)) deallocate(dq)
  allocate (dq(nspecies))
  if (allocated(Q0_TOT)) deallocate(Q0_TOT)
  allocate (Q0_TOT(natoms))
  if (allocated(Q_partial)) deallocate(Q_partial)
  allocate (Q_partial(natoms))
  if (allocated(dq_DP)) deallocate(dq_DP)
  allocate (dq_DP(natoms))

  call initboxes ()

  ! Count the orbitals
  norbitals = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    norbitals = norbitals + num_orb(in1)
  end do
  norbitals_new = norbitals

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
!  Qin = 0.0d0
!  do iatom = 1, natoms
!    in1 = imass(iatom)
!    do issh = 1, nssh(in1)
!      Qin(issh,iatom) = Qneutral(issh,in1)
!    end do
!  end do


  ztot = real(qstate, double)
  nelectron = 0.0d0
  do iatom = 1, natoms
   in1 = imass(iatom)
   do issh = 1, nssh(in1)
    ztot = ztot + Qneutral(issh,in1)
    nelectron(iatom) = nelectron(iatom) + nint(Qneutral(issh,in1))
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

  if (allocated(getmssh)) deallocate(getmssh)
  allocate (getmssh(norbitals))
  if (allocated(getlssh)) deallocate(getlssh)
  allocate (getlssh(norbitals))
  if (allocated(getissh)) deallocate(getissh)
  allocate (getissh(norbitals))
  if (allocated(getiatom)) deallocate(getiatom)
  allocate (getiatom(norbitals))
  if (allocated(get_shell_oforb)) deallocate(get_shell_oforb)
  allocate (get_shell_oforb(norbitals))
  if (allocated(get_orb_ofshell)) deallocate(get_orb_ofshell)
  allocate (get_orb_ofshell(nssh_tot))
  if (allocated(get_issh_ofshell)) deallocate(get_issh_ofshell)
  allocate (get_issh_ofshell(nssh_tot))
  if (allocated(get_l_ofshell)) deallocate(get_l_ofshell)
  allocate (get_l_ofshell(nssh_tot))
  if (allocated(get_iatom_ofshell)) deallocate(get_iatom_ofshell)
  allocate (get_iatom_ofshell(nssh_tot))
 

   alpha = 0 
   imu=1
   print*,'alpha,imu,iatom,issh,lssh'
   do iatom=1,natoms
    do issh=1,nssh(imass(iatom))
      alpha = alpha + 1
      print '(5I6)', alpha, imu, iatom, issh, lssh(issh, imass(iatom))      
      get_orb_ofshell(alpha) = imu
      get_iatom_ofshell(alpha) = iatom
      get_issh_ofshell(alpha) = issh
      get_l_ofshell(alpha) = lssh(issh,imass(iatom))
      imu = imu + ( 2*lssh(issh,imass(iatom))+1 )
     end do
    end do

           
  imu = 0
  alpha = 0 
  do iatom=1,natoms
    do issh=1,nssh(imass(iatom))
      alpha = alpha + 1
      do iorb = 1, 2*lssh(issh,imass(iatom))+1 
        imu=imu+1
        getissh(imu)=issh
        getlssh(imu)=lssh(issh,imass(iatom))
        getiatom(imu)=iatom
        getmssh(imu)=iorb
        get_shell_oforb(imu)=alpha
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
    neigh_max = max(neigh_max, num_neigh)
  end do


 
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

  if (allocated(neigh_b)) deallocate(neigh_b)
  allocate (neigh_b (neigh_max, natoms))
  if (allocated(neigh_j)) deallocate(neigh_j)
  allocate (neigh_j (neigh_max, natoms))
  if (allocated(neigh_comb)) deallocate(neigh_comb)
  allocate (neigh_comb (2, neigh_max**2, natoms))
  if (allocated(neigh_comj)) deallocate(neigh_comj)
  allocate (neigh_comj (2, neigh_max**2, natoms))
  if (allocated(neigh_comm)) deallocate(neigh_comm)
  allocate (neigh_comm (neigh_max**2, natoms))
  if (allocated(neigh_comn)) deallocate(neigh_comn)
  allocate (neigh_comn (natoms))
  if (allocated(neigh_back)) deallocate(neigh_back)
  allocate (neigh_back (natoms, neigh_max))
  if (allocated(neigh_self)) deallocate(neigh_self)
  allocate (neigh_self (natoms))
  if (allocated(nPP_b)) deallocate(nPP_b)
  allocate (nPP_b (neighPP_max, natoms))
  if (allocated(nPP_j)) deallocate(nPP_j)
  allocate (nPP_j (neighPP_max, natoms))
  if (allocated(nPP_map)) deallocate(nPP_map)
  allocate (nPP_map (neighPP_max, natoms))
  if (allocated(nPPn)) deallocate(nPPn)
  allocate (nPPn (natoms))
  if (allocated(nPP_self)) deallocate(nPP_self)
  allocate (nPP_self (natoms))
  if (allocated(nPPx_b)) deallocate(nPPx_b)
  allocate (nPPx_b (neighPP_max, natoms))
  if (allocated(nPPx_j)) deallocate(nPPx_j)
  allocate (nPPx_j (neighPP_max, natoms))
  if (allocated(nPPx_map)) deallocate(nPPx_map)
  allocate (nPPx_map (neighPP_max, natoms))
  if (allocated(nPPx_point)) deallocate(nPPx_point)
  allocate (nPPx_point (neighPP_max, natoms))
  if (allocated(nPPxn)) deallocate(nPPxn)
  allocate (nPPxn (natoms))
  if (allocated(neigh_pair_a1)) deallocate(neigh_pair_a1)
  allocate (neigh_pair_a1 (neigh_max*natoms))
  if (allocated(neigh_pair_a2)) deallocate(neigh_pair_a2)
  allocate (neigh_pair_a2 (neigh_max*natoms))
  if (allocated(neigh_pair_n1)) deallocate(neigh_pair_n1)
  allocate (neigh_pair_n1 (neigh_max*natoms))
  if (allocated(neigh_pair_n2)) deallocate(neigh_pair_n2)
  allocate (neigh_pair_n2 (neigh_max*natoms))
  if (allocated(nPPx_self)) deallocate(nPPx_self)
  allocate (nPPx_self (natoms))
  if (allocated(vxc_1c)) deallocate(vxc_1c)
  allocate (vxc_1c (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(neighn)) deallocate(neighn)
  allocate (neighn (natoms))
  if (allocated(neighn_tot)) deallocate(neighn_tot)
  allocate (neighn_tot (natoms))

! neighPP
  if (allocated(neighPP_b)) deallocate(neighPP_b)
  allocate (neighPP_b (neighPP_max**2, natoms))
  if (allocated(neighPP_j)) deallocate(neighPP_j)
  allocate (neighPP_j (neighPP_max**2, natoms))
  if (allocated(neighPPn)) deallocate(neighPPn)
  allocate (neighPPn (natoms))
! 3. party PPcommon pairs
  if (allocated(neighPP_comb)) deallocate(neighPP_comb)
  allocate (neighPP_comb (2, neighPP_max**2, natoms))
  if (allocated(neighPP_comj)) deallocate(neighPP_comj)
  allocate (neighPP_comj (2, neighPP_max**2, natoms))
  if (allocated(neighPP_comm)) deallocate(neighPP_comm)
  allocate (neighPP_comm (neighPP_max**2, natoms))
  if (allocated(neighPP_comn)) deallocate(neighPP_comn)
  allocate (neighPP_comn (natoms))
  if (allocated(neighPP_self)) deallocate(neighPP_self)
  allocate (neighPP_self (natoms))

! Total neighbor list (mapping together neigh and neighPP part)
  if (allocated(neighj_tot)) deallocate(neighj_tot)
  allocate (neighj_tot (neigh_max+neighPP_max, natoms))
  if (allocated(neighb_tot)) deallocate(neighb_tot)
  allocate (neighb_tot (  neigh_max+neighPP_max, natoms))

  
 
  if (allocated(ioccupy)) deallocate(ioccupy)
  allocate(ioccupy(norbitals))
  if (allocated(ioccupy_k)) deallocate(ioccupy_k)
  allocate(ioccupy_k (norbitals, nkpoints))
  if (allocated(foccupy)) deallocate(foccupy)
  allocate(foccupy (norbitals, nkpoints))
  if (allocated(eigen_k)) deallocate(eigen_k)
  allocate (eigen_k (norbitals, nkpoints))

  
  if (allocated(dewald)) deallocate(dewald)
  allocate (dewald (3, natoms, natoms))
  if (allocated(fewald)) deallocate(fewald)
  allocate (fewald (3, natoms))
  if (allocated(ewald)) deallocate(ewald)
  allocate (ewald (natoms, natoms))


  if (allocated(s_mat)) deallocate(s_mat)
  allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(t_mat)) deallocate(t_mat)
  allocate (t_mat (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(h_mat)) deallocate(h_mat)
  allocate (h_mat (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(sp_mat)) deallocate(sp_mat)
  allocate (sp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(tp_mat)) deallocate(tp_mat)
  allocate (tp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(dipcm)) deallocate(dipcm)
  allocate (dipcm (3, numorb_max, numorb_max))
  if (allocated(dipc)) deallocate(dipc)
  allocate (dipc (3, numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(dip)) deallocate(dip)
  allocate (dip (numorb_max, numorb_max, neigh_max, natoms)) 
  if (allocated(dippcm)) deallocate(dippcm)
  allocate (dippcm (3, 3, numorb_max, numorb_max))
  if (allocated(dippc)) deallocate(dippc)
  allocate (dippc (3, 3, numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(dipp)) deallocate(dipp)
  allocate (dipp (3, numorb_max, numorb_max, neigh_max, natoms)) 
  if (allocated(sm_mat)) deallocate(sm_mat)
  allocate (sm_mat (nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(ewaldqmmm)) deallocate(ewaldqmmm)
  allocate (ewaldqmmm (numorb_max, numorb_max, neigh_max,natoms))
  if (allocated(cape)) deallocate(cape)
  allocate (cape (numorb_max, numorb_max, neigh_max, natoms))

  if (allocated(rho)) deallocate(rho)
  allocate (rho (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(rhoPP)) deallocate(rhoPP)
  allocate (rhoPP (numorb_max, numorb_max, neighPP_max**2, natoms))
  if (allocated(dusr)) deallocate(dusr)
  allocate (dusr (3, natoms))
  

  if (allocated(arho_off)) deallocate(arho_off)
  allocate (arho_off (nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(arhoij_off)) deallocate(arhoij_off)
  allocate (arhoij_off (nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(rho_off)) deallocate(rho_off)
  allocate (rho_off (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(rhoij_off)) deallocate(rhoij_off)
  allocate (rhoij_off (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(arho_on)) deallocate(arho_on)
  allocate (arho_on (nsh_max, nsh_max, natoms))
  if (allocated(arhoi_on)) deallocate(arhoi_on)
  allocate (arhoi_on (nsh_max, nsh_max, natoms))
  if (allocated(rho_on)) deallocate(rho_on)
  allocate (rho_on (numorb_max, numorb_max, natoms))
  if (allocated(rhoi_on)) deallocate(rhoi_on)
  allocate (rhoi_on (numorb_max, numorb_max, natoms))
  if (allocated(vxc_ca)) deallocate(vxc_ca)
  allocate (vxc_ca (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(vxc)) deallocate(vxc)
  allocate (vxc (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(vca)) deallocate(vca)
  allocate (vca (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(ewaldsr)) deallocate(ewaldsr)
  allocate (ewaldsr (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(ewaldlr)) deallocate(ewaldlr)
  allocate (ewaldlr (numorb_max, numorb_max, neigh_max, natoms))

  if (allocated(spm_mat)) deallocate(spm_mat)
  allocate (spm_mat (3, nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(arhop_off)) deallocate(arhop_off)
  allocate (arhop_off (3, nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(arhopij_off)) deallocate(arhopij_off)
  allocate (arhopij_off (3, nsh_max, nsh_max, neigh_max, natoms)) 
  if (allocated(rhop_off)) deallocate(rhop_off)
  allocate (rhop_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
  if (allocated(rhopij_off)) deallocate(rhopij_off)
  allocate (rhopij_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
  if (allocated(arhop_on)) deallocate(arhop_on)
  allocate (arhop_on (3, nsh_max, nsh_max, neigh_max, natoms))
  if (allocated(rhop_on)) deallocate(rhop_on)
  allocate (rhop_on (3, numorb_max, numorb_max, neigh_max, natoms)) 
  ! OLSXC double count corr forces
  if (allocated(dxcdcc)) deallocate(dxcdcc)
  allocate (dxcdcc (3, neigh_max, natoms))

  if (allocated(vna)) deallocate(vna)
  allocate (vna (numorb_max, numorb_max, neigh_max, natoms))
  if (allocated(vnl)) deallocate(vnl)
  allocate (vnl (numorb_max, numorb_max, neighPP_max**2, natoms))
  if (allocated(sVNL)) deallocate(sVNL)
  allocate (sVNL (numorb_max, numorb_max, neighPP_max, natoms))
  if (allocated(spVNL)) deallocate(spVNL)
  allocate (spVNL (3, numorb_max, numorb_max, neighPP_max, natoms))

  call initamat()
  
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

  if (allocated(fana)) deallocate(fana)
  allocate (fana (3, neigh_max, natoms))
  if (allocated(fotna)) deallocate(fotna)
  allocate (fotna (3, neigh_max, natoms)) 
  if (allocated(ft)) deallocate(ft)
  allocate (ft (3, natoms))
  if (allocated(fanl)) deallocate(fanl)
  allocate (fanl (3, neighPP_max, natoms))
  if (allocated(fotnl)) deallocate(fotnl)
  allocate (fotnl (3, neighPP_max, natoms))
  if (allocated(faxc_ca)) deallocate(faxc_ca)
  allocate (faxc_ca (3, neigh_max, natoms)) 
  if (allocated(faxc)) deallocate(faxc)
  allocate (faxc (3, neigh_max, natoms))
  if (allocated(fotxc_ca)) deallocate(fotxc_ca)
  allocate (fotxc_ca (3, neigh_max, natoms))
  if (allocated(fotxc)) deallocate(fotxc)
  allocate (fotxc (3, neigh_max, natoms))
  if (allocated(faca)) deallocate(faca)
  allocate (faca (3, neigh_max, natoms))  
  if (allocated(fotca)) deallocate(fotca)
  allocate (fotca (3, neigh_max, natoms))
  if (allocated(f3naa)) deallocate(f3naa)
  allocate (f3naa (3, natoms))  
  if (allocated(f3nab)) deallocate(f3nab)
  allocate (f3nab (3, natoms))  
  if (allocated(f3nac)) deallocate(f3nac)
  allocate (f3nac (3, natoms))
  if (allocated(f3nla)) deallocate(f3nla)
  allocate (f3nla (3, natoms))  
  if (allocated(f3nlb)) deallocate(f3nlb)
  allocate (f3nlb (3, natoms))  
  if (allocated(f3nlc)) deallocate(f3nlc)
  allocate (f3nlc (3, natoms))
  if (allocated(f3caa)) deallocate(f3caa)
  allocate (f3caa (3, natoms))
  if (allocated(f3cab)) deallocate(f3cab)
  allocate (f3cab (3, natoms))
  if (allocated(f3cac)) deallocate(f3cac)
  allocate (f3cac (3, natoms))
  if (allocated(flrew)) deallocate(flrew)
  allocate (flrew (3, natoms))
  if (allocated(f3xca)) deallocate(f3xca)
  allocate (f3xca (3, natoms))
  if (allocated(f3xcc)) deallocate(f3xcc)
  allocate (f3xcc (3, natoms))
  if (allocated(f3xcb)) deallocate(f3xcb)
  allocate (f3xcb (3, natoms))
  if (allocated(f3xca_ca)) deallocate(f3xca_ca)
  allocate (f3xca_ca (3, natoms))
  if (allocated(f3xcb_ca)) deallocate(f3xcb_ca)
  allocate (f3xcb_ca (3, natoms))
  if (allocated(f3xcc_ca)) deallocate(f3xcc_ca)
  allocate (f3xcc_ca (3, natoms))
  if (allocated(flrew_qmmm)) deallocate(flrew_qmmm)
  allocate (flrew_qmmm (3, natoms))
  if (allocated(fro)) deallocate(fro)
  allocate (fro (3, natoms))
  if (allocated(ftot)) deallocate(ftot)
  allocate (ftot (3, natoms))
  if (allocated(dxcv)) deallocate(dxcv)
  allocate (dxcv (3, natoms))  


  ! kspace variables
  if (allocated(bbnkre)) deallocate(bbnkre)
  allocate (bbnkre (norbitals, norbitals, nkpoints))
  if (icluster .eq. 0 .and. igamma .eq. 0) then
    if (allocated(bbnkim)) deallocate(bbnkim)
    allocate (bbnkim (norbitals, norbitals, nkpoints))
  end if

  if (allocated(g_h)) deallocate(g_h)
  allocate (g_h(numorb_max,numorb_max,nsh_max,natoms,neigh_max,natoms))
  if (allocated(g_xc)) deallocate(g_xc)
  allocate (g_xc(numorb_max,numorb_max,nsh_max,natoms,neigh_max,natoms))
  if (allocated(f_xc)) deallocate(f_xc)
  allocate (f_xc(numorb_max,nsh_max,natoms,neigh_max,natoms))
  if (allocated(exc_aa)) deallocate(exc_aa)
  allocate (exc_aa(numorb_max,natoms,neigh_max,natoms))
  if (allocated(vxc_aa)) deallocate(vxc_aa)
  allocate (vxc_aa(numorb_max,natoms,neigh_max,natoms))

  if (allocated(g_h_shell)) deallocate(g_h_shell)
  allocate (g_h_shell(nsh_max,nsh_max))
  if (allocated(g_xc_shell)) deallocate(g_xc_shell)
  allocate (g_xc_shell(nsh_max,nsh_max))
  if (allocated(f_xc_shell)) deallocate(f_xc_shell)
  allocate (f_xc_shell(nsh_max,nsh_max))
  if (allocated(exc_aa_shell)) deallocate(exc_aa_shell)
  allocate (exc_aa_shell(nsh_max,nsh_max))
  if (allocated(vxc_aa_shell)) deallocate(vxc_aa_shell)
  allocate (vxc_aa_shell(nsh_max,nsh_max))
end subroutine allocate_system
