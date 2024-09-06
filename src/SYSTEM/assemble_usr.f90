subroutine assemble_usr ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system
  use M_fdata
  integer(c_long) iatom
  integer(c_long) ideriv
  integer(c_long) in1, in2
  integer(c_long) index
  integer(c_long) index_coulomb
  integer(c_long) ineigh
  integer(c_long) interaction
  integer(c_long) issh
  integer(c_long) jatom
  integer(c_long) jssh
  integer(c_long) mbeta
  integer(c_long) n1, n2
  real(c_double) distance
  real(c_double) dq1
  real(c_double) dqi
  real(c_double) eklr
  real(c_double) qi, qj
  real(c_double) QQ
  real(c_double) u0tot
  real(c_double) ue0tot
  real(c_double) xforce
  real(c_double) Zi, Zj
  real(c_double), dimension (natoms, neigh_max) :: corksr
  real(c_double), dimension (nsh_max, nsh_max) :: coulomb
  real(c_double), dimension (nsh_max, nsh_max) :: coulombD
  real(c_double), dimension (3) :: dcorksr
  real(c_double), dimension (ME2c_max) :: dslist
  real(c_double), dimension (3) :: eta
  real(c_double), dimension (natoms) :: Q, Q0
  real(c_double), dimension (3) :: r1, r2
  real(c_double), dimension (ME2c_max) :: slist
  real(c_double), dimension (natoms, neigh_max) :: u0
  real(c_double), dimension (natoms) :: uee00
  dusr = 0.0d0
  dxcv = 0.0d0
  u0 = 0.0d0

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
      distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2    + (r2(3) - r1(3))**2)
      index_coulomb = nssh(in1)*nssh(in2)
      interaction = 12
      ideriv = 0
      do index = 1, index_coulomb
        call interpolate_1d (interaction, ideriv, in1, in2, index, iforce, distance, slist(index), dslist(index))
      end do
      n1 = nssh(in1)
      n2 = nssh(in2)

      call recoverC (n1, n2, slist, dslist, coulomb, coulombD)

      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        uee00(iatom) = 0.0d0
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in1)
            uee00(iatom) = uee00(iatom) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        uee00(iatom) = uee00(iatom)*(eq2/2.0d0)
        u0(iatom,ineigh) = 0.0d0
        corksr(iatom,ineigh) = 0.0d0
      else
        u0(iatom,ineigh) = 0.0d0
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            u0(iatom,ineigh) = u0(iatom,ineigh) + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
          end do
        end do
        u0(iatom,ineigh)=(eq2/2.0d0)*(Zi*Zj/distance-u0(iatom,ineigh))
        corksr(iatom,ineigh) = - (eq2/2.0d0)*QQ/distance
        if (iforce .eq. 1) then
          eta(:) = (r2(:) - r1(:))/distance
          xforce = 0.0d0
          do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
              xforce = xforce + Qin(issh,iatom)*Qin(jssh,jatom)*coulombD(issh,jssh)
            end do
          end do
          dusr(:,iatom) = dusr(:,iatom) -  eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
          dusr(:,jatom) = dusr(:,jatom) +  eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
          dcorksr(:) = - eta(:)*(eq2/2.0d0)*QQ/distance**2
          dusr(:,iatom) = dusr(:,iatom) - dcorksr(:)
          dusr(:,jatom) = dusr(:,jatom) + dcorksr(:)

        end if ! forces
      end if ! (iatom .eq. jatom)
    end do
  end do !iatom

  do iatom = 1, natoms
    dusr(:,iatom) = dusr(:,iatom) - (eq2/2.0d0)*fewald(:,iatom)
  end do

  u0tot = 0.0d0
  ue0tot = 0.0d0
  do iatom = 1, natoms
    ue0tot = ue0tot + uee00(iatom)
    do ineigh = 1, neighn(iatom)
      u0tot = u0tot + u0(iatom,ineigh)
      u0tot = u0tot + corksr(iatom,ineigh)
    end do                                     
  end do
  eklr = 0.0d0
  do iatom = 1, natoms
    do jatom = iatom, natoms
     QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
     if (iatom .eq. jatom) then
       eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ
     else
       eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ + (eq2/2.0d0)*ewald(jatom,iatom)*QQ
     end if
    end do
  end do
  u0tot = u0tot - eklr
  uiiuee = u0tot - ue0tot
end subroutine assemble_usr

