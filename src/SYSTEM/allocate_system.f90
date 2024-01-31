subroutine allocate_system ()
  use M_system  

  integer :: iatom, jatom, mbeta, num_neigh, in1,imu,ini2,
  real :: rcutoff_i,distance2,range2

  allocate (ratom (3, natoms))
  allocate (degelec (natoms))
  allocate (imass (natoms))
  allocate (symbol (natoms))
   
  ! allocate_neigh  call allocate_neigh (icluster)
  ! find_neigh_max 
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

        distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
        &         + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2  &
        &         + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2

        ! Add a small displacement to the sum of cutoffs. 
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
        if (distance2 .le. range2) num_neigh = num_neigh + 1
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
        ! Find the distance from (mbeta,jatom) to (0,iatom)
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
        ! Find the distance from (mbeta,jatom) to (0,iatom)
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

    ! Maximum number of neighbors thus far.
    neighPP_max = max(neighPP_max, num_neigh)
  end do


  allocate (neigh_b (neigh_max, natoms))
  allocate (neigh_j (neigh_max, natoms))
  allocate (neighn (natoms))
  allocate (neigh_comb (2, neigh_max**2, natoms))
  allocate (neigh_comj (2, neigh_max**2, natoms))
  if (itheory_xc .eq. 4) allocate (neigh_com_ng (2, neigh_max**2, natoms))
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
  allocate (nPPx_self (natoms))
  allocate (neigh_pair_a1 (neigh_max*natoms))
  allocate (neigh_pair_a2 (neigh_max*natoms))
  allocate (neigh_pair_n1 (neigh_max*natoms))
  allocate (neigh_pair_n2 (neigh_max*natoms))
  if (itheory_xc .eq. 4) allocate (neigh_com_ng (2, neigh_max**2, natoms))
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
  allocate (nPPx_self (natoms))
  allocate (neigh_pair_a1 (neigh_max*natoms))
  allocate (neigh_pair_a2 (neigh_max*natoms))
  allocate (neigh_pair_n1 (neigh_max*natoms))
  allocate (neigh_pair_n2 (neigh_max*natoms))

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
  allocate (neighb_tot (neigh_max+neighPP_max, natoms))
  allocate (neighn_tot (natoms))

! call allocate_f (natoms, neigh_max, neighPP_max, numorb_max, nsh_max, itheory, itheory_xc )
! call allocate_h (natoms, neigh_max, neighPP_max, itheory, itheory_xc,
!call allocate_rho (natoms, neigh_max, neighPP_max, numorb_max,       


!         allocate (arho_off (nsh_max, nsh_max, neigh_max, natoms))
!         allocate (arhoij_off (nsh_max, nsh_max, neigh_max, natoms))
!         allocate (rho_off (numorb_max, numorb_max, neigh_max, natoms))
!         allocate (rhoij_off (numorb_max, numorb_max, neigh_max, natoms))
!         allocate (arho_on (nsh_max, nsh_max, natoms))
!         allocate (arhoi_on (nsh_max, nsh_max, natoms))
!         allocate (rho_on (numorb_max, numorb_max, natoms))
!         allocate (rhoi_on (numorb_max, numorb_max, natoms))

        !call allocate_dos (natoms, iwrtdos, iwrthop)
end subroutine
