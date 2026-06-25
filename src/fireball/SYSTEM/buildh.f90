subroutine buildh ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, vxc_1c, ewaldsr, neigh_j, neighn, vxc, vca, ewaldlr, h_mat, t_mat, vna, ewaldqmmm, &
    & h_mat0, g_h, Qin, Kscf, nssh_tot, get_issh_ofshell, get_iatom_ofshell
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer :: iatom, jatom, imu, inu, ineigh, in1, in2
  integer :: alpha, issh_a, iatom_a
  real(double) :: hmodel, dmax, gcoul, hcoul, dmax_coul
  real(double), dimension(nssh_tot) :: dQ

  h_mat = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          h_mat(imu,inu,ineigh,iatom) = t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + &
            & vxc(imu,inu,ineigh,iatom) + vxc_1c(imu,inu,ineigh,iatom) + vca(imu,inu,ineigh,iatom) + &
            & ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom) + &
            & ewaldqmmm(imu,inu,ineigh,iatom)
        end do
      end do
    end do
  end do

  ! ---- Modelo de cargas estaticas (verificacion no destructiva, todas las proyecciones) ----
  !   H_model = (t + vna)                                 indep. de Q
  !           + vxc(Q) + vxc_1c(Q)                        XC: reconstruido cada paso SCF (no lineal)
  !           + sum_alpha g_h^alpha (Q_alpha - Q^0_alpha) Coulomb/Ewald: congelado en Kscf=1 (lineal, exacto)
  ! En Kscf=1 (dQ=0) guardamos h^0 = h_mat. En Kscf>1 comparamos H_model con el h_mat reensamblado.
  if (Kscf .eq. 1) then
    h_mat0 = h_mat
  else
    do alpha = 1, nssh_tot
      issh_a  = get_issh_ofshell(alpha)
      iatom_a = get_iatom_ofshell(alpha)
      dQ(alpha) = Qin(issh_a,iatom_a) - Qneutral(issh_a,imass(iatom_a))
    end do
    dmax = 0.0d0
    dmax_coul = 0.0d0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            gcoul = dot_product(g_h(:,imu,inu,ineigh,iatom), dQ)
            ! H_model: Coulomb via g_h congelado + XC vivo + (t+vna) + qmmm
            hmodel = t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) &
              & + vxc(imu,inu,ineigh,iatom) + vxc_1c(imu,inu,ineigh,iatom) &
              & + ewaldqmmm(imu,inu,ineigh,iatom) + gcoul
            dmax = max(dmax, abs(hmodel - h_mat(imu,inu,ineigh,iatom)))
            ! parte Coulomb/Ewald aislada: g_h.dQ vs vca+ewaldlr-ewaldsr (debe ser ~0)
            hcoul = vca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom)
            dmax_coul = max(dmax_coul, abs(gcoul - hcoul))
          end do
        end do
      end do
    end do
    print *, 'check H_model vs h_mat  : max|dif| =', dmax,      ' Kscf=', Kscf
    print *, 'check g_h.dQ vs Coulomb : max|dif| =', dmax_coul, ' Kscf=', Kscf
  end if
end subroutine buildh
