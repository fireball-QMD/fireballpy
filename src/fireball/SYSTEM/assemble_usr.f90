subroutine assemble_usr (iauxforce)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, neigh_max, uiiuee, ewald, fewald, neigh_b, neigh_j, neighn, Qin, dq, xl, dusr, dxcv, &
    & Kscf, coulomb_mat, iscf_fast
  use M_fdata, only: nsh_max, ME2c_max, nssh, Qneutral, TWOCENTER_COULOMB, splineint_2c, ind2c
  implicit none
  integer, intent (in) :: iauxforce
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
  real(dp) distance
  real(dp) dq1
  real(dp) dqi
  real(dp) eklr
  real(dp) qi, qj
  real(dp) QQ
  real(dp) u0tot
  real(dp) ue0tot
  real(dp) xforce
  real(dp) Zi, Zj
  real(dp), dimension (natoms, neigh_max) :: corksr
  real(dp), dimension (nsh_max, nsh_max) :: coulomb
  real(dp), dimension (nsh_max, nsh_max) :: coulombD
  real(dp), dimension (3) :: dcorksr
  real(dp), dimension (ME2c_max) :: dslist
  real(dp), dimension (3) :: eta
  real(dp), dimension (natoms) :: Q, Q0
  real(dp), dimension (3) :: r1, r2
  real(dp), dimension (ME2c_max) :: slist
  real(dp), dimension (natoms, neigh_max) :: u0
  real(dp), dimension (natoms) :: uee00
  dusr = 0.0_dp
  dxcv = 0.0_dp
  u0 = 0.0_dp

  do iatom = 1, natoms
    Q(iatom) = 0.0_dp
    Q0(iatom) = 0.0_dp
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
      n1 = nssh(in1)
      n2 = nssh(in2)
      ! Speed up SCF loop: las integrales de Coulomb son Q-independientes; en Kscf>1
      ! (sin fuerzas) se reutiliza la tabla congelada en Kscf=1 en vez de reinterpolar.
      if (Kscf .gt. 1 .and. iauxforce .eq. 0 .and. iscf_fast .eq. 1) then
        coulomb(1:n1,1:n2) = coulomb_mat(1:n1,1:n2,ineigh,iatom)
      else
        index_coulomb = nssh(in1)*nssh(in2)
        interaction = TWOCENTER_COULOMB
        ideriv = 0
        do index = 1, index_coulomb
          slist(index) = splineint_2c(index, ind2c(interaction, ideriv), in1, in2)%f(distance)
          if (iauxforce == 1) dslist(index) = splineint_2c(index, ind2c(interaction, ideriv), in1, in2)%f(distance, order=1)
        end do
        call recoverC (n1, n2, slist, dslist, coulomb, coulombD)
        if (Kscf .eq. 1) coulomb_mat(1:n1,1:n2,ineigh,iatom) = coulomb(1:n1,1:n2)
      end if

      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        uee00(iatom) = 0.0_dp
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in1)
            uee00(iatom) = uee00(iatom) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        uee00(iatom) = uee00(iatom)*(eq2/2.0_dp)
        u0(iatom,ineigh) = 0.0_dp
        corksr(iatom,ineigh) = 0.0_dp
      else
        u0(iatom,ineigh) = 0.0_dp
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            u0(iatom,ineigh) = u0(iatom,ineigh) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        u0(iatom,ineigh)=(eq2/2.0_dp)*(Zi*Zj/distance-u0(iatom,ineigh))
        corksr(iatom,ineigh) = - (eq2/2.0_dp)*QQ/distance
        if (iauxforce .eq. 1) then
          eta(:) = (r2(:) - r1(:))/distance
          xforce = 0.0_dp
          do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
              xforce = xforce + Qin(issh,iatom)*Qin(jssh,jatom)*coulombD(issh,jssh)
            end do
          end do
          dusr(:,iatom) = dusr(:,iatom) -  eta(:)*(eq2/2.0_dp)*(Zi*Zj/distance**2 + xforce)
          dusr(:,jatom) = dusr(:,jatom) +  eta(:)*(eq2/2.0_dp)*(Zi*Zj/distance**2 + xforce)
          dcorksr(:) = - eta(:)*(eq2/2.0_dp)*QQ/distance**2
          dusr(:,iatom) = dusr(:,iatom) - dcorksr(:)
          dusr(:,jatom) = dusr(:,jatom) + dcorksr(:)

        end if ! forces
      end if ! (iatom .eq. jatom)
    end do
  end do !iatom

  do iatom = 1, natoms
    dusr(:,iatom) = dusr(:,iatom) - (eq2/2.0_dp)*fewald(:,iatom)
  end do

  u0tot = 0.0_dp
  ue0tot = 0.0_dp
  do iatom = 1, natoms
    ue0tot = ue0tot + uee00(iatom)
    do ineigh = 1, neighn(iatom)
      u0tot = u0tot + u0(iatom,ineigh)
      u0tot = u0tot + corksr(iatom,ineigh)
    end do                                     
  end do
  eklr = 0.0_dp
  do iatom = 1, natoms
    do jatom = iatom, natoms
     QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
     if (iatom .eq. jatom) then
       eklr = eklr + (eq2/2.0_dp)*ewald(iatom,jatom)*QQ
     else
       eklr = eklr + (eq2/2.0_dp)*ewald(iatom,jatom)*QQ + (eq2/2.0_dp)*ewald(jatom,iatom)*QQ
     end if
    end do
  end do
  u0tot = u0tot - eklr
  uiiuee = u0tot - ue0tot
end subroutine assemble_usr
