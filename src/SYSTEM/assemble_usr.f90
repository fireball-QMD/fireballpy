subroutine assemble_usr ()
  use M_system
  use M_fdata
  use M_constants
  implicit none
  integer iatom
  integer ideriv
  integer in1, in2
  integer index
  integer index_coulomb
  integer ineigh
  integer interaction
  integer issh
  integer jatom
  integer jssh
  integer mbeta
  integer n1, n2
  integer non2c
  real distance
  real dq1, dq2
  real dqi, dqj
  real dxc
  real dxc00, dxc0P, dxc0M, dxcP0, dxcM0
  real eklr
  real qi, qj
  real QQ
  real u0tot
  real ue0tot
  real xc
  real xc00, xc0P, xc0M, xcP0, xcM0
  real xforce
  real Zi, Zj
 
  real, dimension (natoms, neigh_max) :: corksr
  real, dimension (nsh_max, nsh_max) :: coulomb
  real, dimension (nsh_max, nsh_max) :: coulombD
  real, dimension (3) :: dcorksr
  real, dimension (ME2c_max) :: dslist
  real, dimension (3) :: eta
  real, dimension (natoms) :: Q, Q0
  real, dimension (3) :: r1, r2
  real, dimension (ME2c_max) :: slist
  real, dimension (natoms, neigh_max) :: u0
  real, dimension (natoms) :: uee00
  
  dxcv = 0.0d0
  dusr = 0.0d0
  u0 = 0.0d0
  uxcdcc = 0.0d0
 
  do iatom = 1, natoms
    Q(iatom) = 0.0d0
    Q0(iatom) = 0.0d0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
     Q(iatom) = Q(iatom) + Qin(issh,iatom)
     Q0(iatom) = Q0(iatom) + Qneutral(issh,in1)
    end do
  end do
   
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    qi = Q(iatom)
    Zi = Q0(iatom)
    dqi = Q(iatom) - Q0(iatom)
    dq1 = dq(in1)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = xl(:,mbeta) + ratom(:,jatom)
      in2 = imass(jatom)
      qj = Q(jatom)
      Zj = Q0(jatom)
      QQ = Zi*Zj - qi*qj
      distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2 + (r2(3) - r1(3))**2)
 
      ! GET COULOMB INTERACTIONS 
      ! Now find the three coulomb integrals need to evaluate the neutral
      ! atom/neutral atom hartree interaction.
      ! Loop over all the non-zero integrals for this interaction:
      index_coulomb = nssh(in1)*nssh(in2)
      interaction = 12
      ideriv = 0
      do index = 1, index_coulomb
        call interpolate_1d (interaction, ideriv, in1, in2, index, distance, slist(index), dslist(index))
      end do
 
      ! We have the data, it is stored in the following way: v(1,1), v(1,2),
      ! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
      ! ordered according to the index_coulomb. Restore this to true 2x2 format:
      n1 = nssh(in1)
      n2 = nssh(in2)
      call recoverC (n1, n2, slist, dslist, coulomb, coulombD)
         
      ! Actually, now we calculate not only the neutral atom contribution,
      ! but also the short-ranged contribution due to the transfer of charge
      ! between the atoms:
      ! (Eii - Eee)neut - SUM(short range)(n(i) + dn(i))*dn(j)*J(i,j),
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
        ! SPECIAL CASE: SELF-INTERACTION
        uee00(iatom) = 0.0d0
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in1)
            uee00(iatom) = uee00(iatom) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        ! put the half and the units in:
        uee00(iatom) = uee00(iatom)*(eq2/2.0d0)
        u0(iatom,ineigh) = 0.0d0
        corksr(iatom,ineigh) = 0.0d0
      else 
        ! BONAFIDE TWO ATOM CASE
        ! Compute u0
        u0(iatom,ineigh) = 0.0d0
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            u0(iatom,ineigh) = u0(iatom,ineigh) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        u0(iatom,ineigh)=(eq2/2.0d0)*(Zi*Zj/distance-u0(iatom,ineigh)) 

        ! This is a correction for the extention of the long-ranged sum to all the
        ! atoms in the system, which we do in order to use the Ewald summation
        ! technique. This energy is to be subtracted from the total, since it is
        ! correctly calculated by integrals, but overcounted by the ewald
        ! contribution. We include the minus HERE, which means we subtract it from
        ! etot.
        corksr(iatom,ineigh) = - (eq2/2.0d0)*QQ/distance
 
        !************
        !    FORCES
        ! ************
        if (iforce .eq. 1) then
          eta(:) = (r2(:) - r1(:))/distance
          ! Put in the forces due to the charge transfer. This 'sumit' has the sign of
          ! d/rd1, and is NOT force-like. We put in force-like character later in dusr.
          xforce = 0.0d0
          if (itheory .eq. 1) then
            do issh = 1, nssh(in1)
              do jssh = 1, nssh(in2)
                xforce = xforce + Qin(issh,iatom)*Qin(jssh,jatom)*coulombD(issh,jssh)
              end do
            end do
            dusr(:,iatom) = dusr(:,iatom) - eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
            dusr(:,jatom) = dusr(:,jatom) + eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
            ! Now we add the corksr correction. Both of these are d/dr1
            ! derivatives and are NOT force-like.
            dcorksr(:) = - eta(:)*(eq2/2.0d0)*QQ/distance**2
            dusr(:,iatom) = dusr(:,iatom) - dcorksr(:)
            dusr(:,jatom) = dusr(:,jatom) + dcorksr(:)
          end if
        end if  ! end if (forces)
      end if  ! end if (iatom .eq. jatom)
    end do  ! End of loop over neighbors
  end do ! End of loop over iatom

  ! Subtract the forces for the ewald interaction
  ! The variable fewald is already force-like.
  if (itheory .eq. 1) then 
    do iatom = 1, natoms
      dusr(:,iatom) = dusr(:,iatom) - (eq2/2.0d0)*fewald(:,iatom)
    end do
  end if
  
  ! ***************************************************************************
  ! Compute the total cell value of uii-uee; uii-uee = sum u0(i,m) - sum uee00(i)
  ! Add long-range ewald interactions.
  ! ***************************************************************************
  u0tot = 0.0d0
  ue0tot = 0.0d0
  do iatom = 1, natoms
    ue0tot = ue0tot + uee00(iatom)
    do ineigh = 1, neighn(iatom)
      u0tot = u0tot + u0(iatom,ineigh)
      if (itheory .eq. 1) u0tot = u0tot + corksr(iatom,ineigh)
    end do
  end do
 
  ! Notice that we add the corksr correction to etot.
  ! This is because of the subtraction minus discussed above.
  eklr = 0.0d0
  if (itheory .eq. 1) then
    do iatom = 1, natoms
      do jatom = iatom, natoms
        ! Calculate q(iatom)*q(jatom) - q0(iatom)*q0(jatom) = QQ
        QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
        if (iatom .eq. jatom) then
          eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ
        else
          eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ + (eq2/2.0d0)*ewald(jatom,iatom)*QQ
        end if
      end do
    end do
    u0tot = u0tot - eklr
  end if
  uiiuee = u0tot - ue0tot
  if (V_intra_dip .eq. 1) then
    !write(*,*) 'Ankais dc_v_intra_dip = ', dc_v_intra_dip_1c
    uiiuee = uiiuee + dc_v_intra_dip_1c
    !Double-counting V_intra_dip_1c here??? Double counting
  end if !end if V_intra_dip .eq. 1

  return
end subroutine assemble_usr
