subroutine buildh ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, vxc_1c, ewaldsr, neigh_j, neighn, vxc, vca, ewaldlr, h_mat, t_mat, vna, ewaldqmmm, &
    & h_mat0, g_h, Qin, Kscf, nssh_tot, get_issh_ofshell, get_iatom_ofshell
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer :: iatom, jatom, imu, inu, ineigh, in1, in2
  integer :: alpha, issh_a, iatom_a
  real(double), dimension(nssh_tot) :: dQ

  ! ---- Cargas estaticas: H_{mu nu} = h^0_{mu nu} + sum_alpha g_h^alpha_{mu nu} (Q_alpha - Q^0_alpha) ----
  ! Kscf=1 (dQ=0): reensamblado completo -> h^0 (en este paso vca=ewaldlr=ewaldsr=0).
  ! Kscf>1       : Coulomb/Ewald via g_h congelado (dot_product), XC reconstruido cada paso.
  if (Kscf .eq. 1) then
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
    h_mat0 = h_mat
  else
    do alpha = 1, nssh_tot
      issh_a  = get_issh_ofshell(alpha)
      iatom_a = get_iatom_ofshell(alpha)
      dQ(alpha) = Qin(issh_a,iatom_a) - Qneutral(issh_a,imass(iatom_a))
    end do
    h_mat = 0.0d0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            h_mat(imu,inu,ineigh,iatom) = t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom) + &
              & vxc(imu,inu,ineigh,iatom) + vxc_1c(imu,inu,ineigh,iatom) + ewaldqmmm(imu,inu,ineigh,iatom) + &
              & dot_product(g_h(:,imu,inu,ineigh,iatom), dQ)
          end do
        end do
      end do
    end do
  end if
end subroutine buildh
