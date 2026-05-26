! This routine calculates the average densities with charge transfer.
subroutine average_ca_rho()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_fdata, only: nssh, num_orb, nsh_max
  use M_system, only: iforce, xc_overtol, natoms, ratom, imass, neigh_max, Kscf, &
  &                   neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, &
  &                   neigh_comm, neigh_comn, neigh_back, numorb_max, Qin, rho_off, &
  &                   rhoij_off, rho_on, arho_on, rhoi_on, &
  &                   arhoi_on, arhop_on, rhop_on, arhoij_off, arho_off, arhopij_off, &
  &                   arhop_off, rhop_off, rhopij_off, xl
  implicit none
  integer :: iatom, ibeta, imu, in1, in2, indna, ineigh, interaction, interaction0, &
  &          inu, isorp, ialp, issh, jatom, jneigh, jbeta, jssh, mbeta, mneigh, interaction3c
  real(dp) :: cost, x, y
  real(dp) :: r1(3), r2(3), r21(3), rhat(3), rnabc(3), rna(3), sighat(3), eps(3, 3), &
  &           deps(3, 3, 3), rhomx(numorb_max, numorb_max), rhompx(3, numorb_max, numorb_max), &
  &           rhomm(nsh_max, nsh_max), rhompm(3, nsh_max, nsh_max), &
  &           rhom_2c(nsh_max, nsh_max, neigh_max, natoms), rhomi_2c(nsh_max, nsh_max, natoms), &
  &           rhomp_2c(3, nsh_max, nsh_max, neigh_max, natoms)

  !   -----  ON SITE PART  ------
  rho_on = 0.0_dp
  rhoi_on = 0.0_dp
  rhop_on = 0.0_dp
  rhom_2c = 0.0_dp
  rhomi_2c = 0.0_dp
  rhomp_2c = 0.0_dp
  interaction = 17
  interaction0 = 22
  do iatom = 1, natoms
    r1(:) = ratom(:, iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh, iatom)
      jatom = neigh_j(ineigh, iatom)
      in2 = imass(jatom)
      r2 = ratom(:, jatom) + xl(:, mbeta)
      r21 = r2 - r1
      y = norm2(r21)
      if (y < 1.0e-05_dp) then
        sighat(1) = 0.0_dp
        sighat(2) = 0.0_dp
        sighat(3) = 1.0_dp
      else
        sighat = r21/y
      end if
      call epsilon(r2, sighat, eps)
      call deps2cent(r1, r2, eps, deps)

      if (iatom == jatom .and. mbeta == 0) then
        do isorp = 1, nssh(in2)
          call doscentros(interaction, isorp, iforce, in1, in2, in1, y, eps, deps, rhomx, rhompx)
          call doscentrosS(interaction0, isorp, iforce, in1, in2, in1, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
              rhoi_on(imu, inu, iatom) = rhoi_on(imu, inu, iatom) + &
              &                                rhomx(imu, inu)*Qin(isorp, jatom)
              rho_on(imu, inu, iatom) = rho_on(imu, inu, iatom) + &
              &                               rhomx(imu, inu)*Qin(isorp, jatom)
            end do
          end do
          do jssh = 1, nssh(in1)
            do issh = 1, nssh(in1)
              rhomi_2c(issh, jssh, iatom) = rhomi_2c(issh, jssh, iatom) + &
              &                                    rhomm(issh, jssh)*Qin(isorp, jatom)
              rhom_2c(issh, jssh, ineigh, iatom) = rhom_2c(issh, jssh, ineigh, iatom) + &
              &                                    rhomm(issh, jssh)*Qin(isorp, jatom)
            end do   ! issh
          end do   ! jssh
        end do   ! isorp
      else
        do isorp = 1, nssh(in2)
          call doscentros(interaction, isorp, iforce, in1, in2, in1, y, eps, deps, rhomx, rhompx)
          call doscentrosS(interaction0, isorp, iforce, in1, in2, in1, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
              rhop_on(:, imu, inu, ineigh, iatom) = rhop_on(:, imu, inu, ineigh, iatom) + &
              &                                 rhompx(:, imu, inu)*Qin(isorp, jatom)
              rho_on(imu, inu, iatom) = rho_on(imu, inu, iatom) + &
              &                               rhomx(imu, inu)*Qin(isorp, jatom)
            end do
          end do
          do jssh = 1, nssh(in1)
            do issh = 1, nssh(in1)
              rhom_2c(issh, jssh, ineigh, iatom) = rhom_2c(issh, jssh, ineigh, iatom) + &
              &                                    rhomm(issh, jssh)*Qin(isorp, jatom)
              rhomp_2c(:, issh, jssh, ineigh, iatom) = rhomp_2c(:, issh, jssh, ineigh, iatom) + &
              &                                        rhompm(:, issh, jssh)*Qin(isorp, jatom)
            end do   ! issh
          end do   ! jssh
        end do   ! isorp
      end if
    end do   ! ineigh
  end do

  ! AVERAGE DENSITY APPROXIMATION (not yet 3c part)
  arho_on = 0.0_dp
  arhoi_on = 0.0_dp
  arhop_on = 0.0_dp
  arhoij_off = 0.0_dp
  arhopij_off = 0.0_dp
  arho_off = 0.0_dp
  arhop_off = 0.0_dp
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh, iatom)
      jatom = neigh_j(ineigh, iatom)
      in2 = imass(jatom)
      do jssh = 1, nssh(in1)
        do issh = 1, nssh(in1)
          arho_on(issh, jssh, iatom) = arho_on(issh, jssh, iatom) + &
            &                          0.5_dp*(rhom_2c(issh, issh, ineigh, iatom) + rhom_2c(jssh, jssh, ineigh, iatom))
          arhop_on(:, issh, jssh, ineigh, iatom) = 0.5_dp*(rhomp_2c(:, issh, issh, ineigh, iatom) + &
            &                                              rhomp_2c(:, jssh, jssh, ineigh, iatom))
        end do
      end do
      if (iatom == jatom .and. mbeta == 0) then
        do jssh = 1, nssh(in1)
          do issh = 1, nssh(in1)
            arhoi_on(issh, jssh, iatom) = 0.5_dp*(rhomi_2c(issh, issh, iatom) + &
              &                                   rhomi_2c(jssh, jssh, iatom))
          end do
        end do
      else
        jneigh = neigh_back(iatom, ineigh)
        do jssh = 1, nssh(in2)
          do issh = 1, nssh(in1)
            arhoij_off(issh, jssh, ineigh, iatom) = 0.5_dp*(rhomi_2c(issh, issh, iatom) + &
              &                                             rhomi_2c(jssh, jssh, jatom) + &
              &                                             rhom_2c(issh, issh, ineigh, iatom) + &
              &                                             rhom_2c(jssh, jssh, jneigh, jatom))
            arhopij_off(:, issh, jssh, ineigh, iatom) = 0.5_dp*(rhomp_2c(:, issh, issh, ineigh, iatom) + &
              &                                                 rhomp_2c(:, jssh, jssh, jneigh, jatom))
          end do
        end do
      end if
    end do   ! do ineigh
  end do    ! do iatom
  arho_off = arhoij_off
  arhop_off = arhopij_off

  ! THREE CENTER PART
  rho_off = 0.0_dp
  rhop_off = 0.0_dp
  interaction3c = 3
  do ialp = 1, natoms
    rna = ratom(:, ialp)
    indna = imass(ialp)
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh, ialp)
      if (mneigh == 0) cycle
      iatom = neigh_comj(1, ineigh, ialp)
      ibeta = neigh_comb(1, ineigh, ialp)
      r1 = ratom(:, iatom) + xl(:, ibeta)
      in1 = imass(iatom)
      jatom = neigh_comj(2, ineigh, ialp)
      jbeta = neigh_comb(2, ineigh, ialp)
      r2 = ratom(:, jatom) + xl(:, jbeta)
      in2 = imass(jatom)
      jneigh = neigh_back(iatom, mneigh)
      r21 = r2 - r1
      y = norm2(r21)
      if (y < 1.0e-05_dp) then
        sighat(1) = 0.0_dp
        sighat(2) = 0.0_dp
        sighat(3) = 1.0_dp
      else
        sighat(:) = r21(:)/y
      end if
      rnabc = rna - (r1 + 0.5_dp*r21)
      x = norm2(rnabc)
      if (x < 1.0e-05_dp) then
        rhat(1) = 0.0_dp
        rhat(2) = 0.0_dp
        rhat(3) = 0.0_dp
      else
        rhat = rnabc/x
      end if
      cost = dot_product(sighat, rhat)
      call epsilon(rhat, sighat, eps)
      do isorp = 1, nssh(indna)
        call trescentros(interaction3c, isorp, in1, in2, indna, x, y, cost, eps, rhomx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rho_off(imu, inu, mneigh, iatom) = rho_off(imu, inu, mneigh, iatom) + &
            &                               rhomx(imu, inu)*Qin(isorp, ialp)
            rho_off(inu, imu, jneigh, jatom) = rho_off(imu, inu, mneigh, iatom)
          end do ! imu
        end do ! inu
      end do ! isorp
      do jssh = 1, nssh(in2)
        do issh = 1, nssh(in1)
          arho_off(issh, jssh, mneigh, iatom) = arho_off(issh, jssh, mneigh, iatom) + &
            &                                   0.5_dp*(rhom_2c(issh, issh, mneigh, iatom) + &
            &                                           rhom_2c(jssh, jssh, jneigh, jatom))
          arhop_off(:, issh, jssh, mneigh, iatom) = arhop_off(:, issh, jssh, mneigh, iatom) + &
            &                                       0.5_dp*rhomp_2c(:, issh, issh, mneigh, iatom)
        end do
      end do
    end do ! ineigh
  end do ! ialp

  !   -----  OFF SITE PART  ------
  ! We assemble off-site density matrices
  ! <mu i| (rho_i+rho_j) | nu j>  ......  rhoij_off; arhoij_off
  ! <mu i| rho_k | nu j>   ......  rho_off; arho_off
  rhoij_off = 0.0_dp
  rhopij_off = 0.0_dp
  do iatom = 1, natoms
    r1(:) = ratom(:, iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh, iatom)
      jatom = neigh_j(ineigh, iatom)
      r2(:) = ratom(:, jatom) + xl(:, mbeta)
      in2 = imass(jatom)
      if (iatom == jatom .and. mbeta == 0) cycle
      r21(:) = r2(:) - r1(:)
      y = norm2(r21)
      if (y < 1.0e-05_dp) then
        sighat(1) = 0.0_dp
        sighat(2) = 0.0_dp
        sighat(3) = 1.0_dp
      else
        sighat(:) = r21(:)/y
      end if
      call epsilon(r2, sighat, eps)
      call deps2cent(r1, r2, eps, deps)
      interaction = 15
      do isorp = 1, nssh(in1)
        call doscentros(interaction, isorp, iforce, in1, in1, in2, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rhoij_off(imu, inu, ineigh, iatom) = rhoij_off(imu, inu, ineigh, iatom) + &
              &                                  rhomx(imu, inu)*Qin(isorp, iatom)
            rhopij_off(:, imu, inu, ineigh, iatom) = rhopij_off(:, imu, inu, ineigh, iatom) + &
              &                                      rhompx(:, imu, inu)*Qin(isorp, iatom)
          end do ! imu
        end do ! inu
      end do ! isorp

      interaction = 16
      do isorp = 1, nssh(in2)
        call doscentros(interaction, isorp, iforce, in1, in2, in2, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rhoij_off(imu, inu, ineigh, iatom) = rhoij_off(imu, inu, ineigh, iatom) + &
              &                                  rhomx(imu, inu)*Qin(isorp, jatom)
            rhopij_off(:, imu, inu, ineigh, iatom) = rhopij_off(:, imu, inu, ineigh, iatom) + &
              &                                      rhompx(:, imu, inu)*Qin(isorp, jatom)
          end do ! imu
        end do ! inu
      end do ! isorp
    end do ! ineigh
  end do ! iatom
  rho_off = rho_off + rhoij_off
  rhop_off = rhop_off + rhopij_off
end subroutine average_ca_rho
