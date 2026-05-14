! This routine calculates the average densities with charge transfer.
subroutine average_ca_rho ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nssh, num_orb, nsh_max
  use M_system, only: iforce, xc_overtol, natoms, ratom, imass, neigh_max, Kscf, &
  &                   neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, &
  &                   neigh_comm, neigh_comn, neigh_back, numorb_max, Qin, rho_off, &
  &                   rhoij_off, rho_on, arho_on, rhoi_on, &
  &                   arhoi_on, arhop_on, rhop_on, arhoij_off, arho_off, arhopij_off, &
  &                   arhop_off, rhop_off, rhopij_off, xl
  implicit none
  integer :: iatom, ibeta, imu, in1, in2, indna, ineigh, interaction, interaction0, &
  &          inu, isorp, ialp, issh, jatom, jneigh, jbeta, jssh, mbeta, mneigh
  real(double) :: cost, x, y
  real(double) :: r1(3), r2(3), r21(3), rhat(3), rnabc(3), rna(3), sighat(3), eps(3,3), &
  &               deps(3,3,3), rhomx(numorb_max, numorb_max), rhompx(3, numorb_max, numorb_max), &
  &               rhomm(nsh_max, nsh_max), rhompm(3, nsh_max, nsh_max), &
  &               rhom_2c(nsh_max, nsh_max, natoms), rhomi_2c(nsh_max, nsh_max, natoms), &
  &               rhomp_2c(3, nsh_max, nsh_max, natoms)

  !   -----  ON SITE PART  ------
  rho_on      = 0.0d0
  rhoi_on     = 0.0d0
  rhop_on     = 0.0d0
  rhom_2c     = 0.0d0
  rhomi_2c    = 0.0d0
  rhomp_2c    = 0.0d0
  arho_on     = 0.0d0
  arhoi_on    = 0.0d0
  arhop_on    = 0.0d0
  arhoij_off  = 0.0d0
  arhopij_off = 0.0d0
  arho_off    = 0.0d0
  arhop_off   = 0.0d0
  do iatom = 1, natoms
    r1(:) = ratom(:, iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      r2 = ratom(:,jatom) + xl(:,mbeta)
      r21 = r2 - r1
      y = norm2(r21)
      if (y < 1.0d-05) then
        sighat(1) = 0.0d0
        sighat(2) = 0.0d0
        sighat(3) = 1.0d0
      else
        sighat = r21/y
      end if
      call epsilon (r2, sighat, eps)
      call deps2cent (r1, r2, eps, deps)

      ! CALL DOSCENTROS AND GET VXC FOR ATM CASE - AVERAGE DENSITY APPROXIMATION
      interaction = 17
      interaction0 = 22
      do isorp = 1, nssh(in2)
        call doscentros (interaction, isorp, iforce, in1, in2, in1, y, eps, deps, rhomx, rhompx)
        call doscentrosS (interaction0, isorp, iforce, in1, in2, in1, y, eps, rhomm, rhompm)
        if (iatom == jatom .and. mbeta == 0) then
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
              rhoi_on(imu,inu,iatom) = rhoi_on(imu,inu,iatom) + &
              &                                rhomx(imu,inu)*Qin(isorp,jatom)
              rho_on(imu,inu,iatom)  = rho_on(imu,inu,iatom)  + &
              &                               rhomx(imu,inu)*Qin(isorp,jatom)
            end do
          end do
          do issh = 1, nssh(in1)
            do jssh = 1, nssh(in1)
              rhomi_2c(issh,jssh,iatom) = rhomi_2c(issh,jssh,iatom) + &
              &                           rhomm(issh,jssh)*Qin(isorp,jatom)
              rhom_2c(issh,jssh,iatom)  = rhom_2c(issh,jssh,iatom)  + &
              &                           rhomm(issh,jssh)*Qin(isorp,jatom)
              arhoi_on(issh,jssh,iatom) = arhoi_on(issh,jssh,iatom) + &
              &                         0.5d0*(rhomm(issh,issh) + rhomm(jssh,jssh))*Qin(isorp,jatom)
            end do   ! issh
          end do   ! jssh
        else
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
              rho_on(imu,inu,iatom)           = rho_on(imu,inu,iatom)           + &
              &                                 rhomx(imu,inu)*Qin(isorp,jatom)
              rhop_on(:,imu,inu,ineigh,iatom) = rhop_on(:,imu,inu,ineigh,iatom) + &
              &                                 rhompx(:,imu,inu)*Qin(isorp,jatom)
            end do
          end do
          do issh = 1, nssh(in1)
            do jssh = 1, nssh(in1)
              rhom_2c(issh,jssh,iatom)    = rhom_2c(issh,jssh,iatom)    + &
              &                             rhomm(issh,jssh)*Qin(isorp,jatom)
              rhomp_2c(:,issh,jssh,iatom) = rhomp_2c(:,issh,jssh,iatom) + &
              &                             rhompm(:,issh,jssh)*Qin(isorp,jatom)
              arhop_on(:,issh,jssh,ineigh,iatom) = arhop_on(:,issh,jssh,ineigh,iatom) + &
              &                   0.5d0*(rhompm(:,issh,issh) + rhompm(:,jssh,jssh))*Qin(isorp,jatom)
            end do   ! jssh
          end do   ! issh
        end if
      end do   ! isorp
    end do   ! ineigh
  end do

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh, iatom)
      jatom = neigh_j(ineigh, iatom)
      in2 = imass(jatom)
      if (iatom == jatom .and. mbeta == 0) then
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in1)
            arho_on(issh,jssh,iatom)  = 0.5d0*(rhom_2c(issh,issh,iatom)  + &
            &                                  rhom_2c(jssh,jssh,iatom))
          end do
        end do
      else
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            arhoij_off(issh,jssh,ineigh,iatom)  = 0.5d0*(rhomi_2c(issh,issh,iatom)   + &
            &                                            rhomi_2c(jssh,jssh,jatom))
            arho_off(issh,jssh,ineigh,iatom)    = 0.5d0*(rhom_2c(issh,issh,iatom)    + &
            &                                            rhom_2c(jssh,jssh,jatom))
            arhop_off(:,issh,jssh,ineigh,iatom) = 0.5d0*(rhomp_2c(:,issh,issh,iatom) + &
            &                                            rhomp_2c(:,jssh,jssh,jatom))
          end do
        end do
      end if
    end do   ! do ineigh
  end do    ! do iatom

  !   -----  OFF SITE PART  ------
  ! We assemble off-site density matrices
  ! <mu i| (rho_i+rho_j) | nu j>  ......  rhoij_off; arhoij_off
  ! <mu i| rho_k | nu j>   ......  rho_off; arho_off

  ! two center
  rho_off    = 0.0d0
  rhop_off   = 0.0d0
  rhoij_off  = 0.0d0
  rhopij_off = 0.0d0
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      if (iatom == jatom .and. mbeta == 0) cycle
      r21(:) = r2(:) - r1(:)
      y = norm2(r21)
      if (y < 1.0d-05) then
        sighat(1) = 0.0d0
        sighat(2) = 0.0d0
        sighat(3) = 1.0d0
      else
        sighat(:) = r21(:)/y
      end if
      call epsilon (r2, sighat, eps)
      call deps2cent (r1, r2, eps, deps)
      interaction = 15
      do isorp = 1, nssh(in1)
        call doscentros (interaction, isorp, iforce, in1, in1, in2, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rhoij_off(imu,inu,ineigh,iatom)    = rhoij_off(imu,inu,ineigh,iatom)    + &
            &                                    rhomx(imu,inu)*Qin(isorp,iatom)
            rhopij_off(:,imu,inu,ineigh,iatom) = rhopij_off(:,imu,inu,ineigh,iatom) + &
            &                                    rhompx(:,imu,inu)*Qin(isorp,iatom)
          end do ! imu
        end do ! inu
      end do ! isorp

      interaction = 16
      do isorp = 1, nssh(in2)
        call doscentros (interaction, isorp, iforce, in1, in2, in2, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rhoij_off(imu,inu,ineigh,iatom)    = rhoij_off(imu,inu,ineigh,iatom)    + &
            &                                    rhomx(imu,inu)*Qin(isorp,jatom)
            rhopij_off(:,imu,inu,ineigh,iatom) = rhopij_off(:,imu,inu,ineigh,iatom) + &
            &                                    rhompx(:,imu,inu)*Qin(isorp,jatom)
          end do ! imu
        end do ! inu
      end do ! isorp
    end do ! ineigh
  end do ! iatom
  rho_off = rhoij_off
  rhop_off = rhopij_off

  ! three center
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh,ialp)
      if (mneigh == 0) cycle
      iatom = neigh_comj(1,ineigh,ialp)
      ibeta = neigh_comb(1,ineigh,ialp)
      r1(:) = ratom(:,iatom) + xl(:,ibeta)
      in1 = imass(iatom)
      jatom = neigh_comj(2,ineigh,ialp)
      jbeta = neigh_comb(2,ineigh,ialp)
      r2(:) = ratom(:,jatom) + xl(:,jbeta)
      in2 = imass(jatom)
      jneigh = neigh_back(iatom,mneigh)
      r21(:) = r2(:) - r1(:)
      y = norm2(r21)
      if (y < 1.0d-05) then
        sighat(1) = 0.0d0
        sighat(2) = 0.0d0
        sighat(3) = 1.0d0
      else
        sighat(:) = r21(:)/y
      end if
      rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
      x = norm2(rnabc)
      if (x < 1.0d-05) then
        rhat(1) = 0.0d0
        rhat(2) = 0.0d0
        rhat(3) = 0.0d0
      else
        rhat(:) = rnabc(:)/x
      end if
      cost = dot_product(sighat, rhat)
      call epsilon (rhat, sighat, eps)
      interaction = 3
      do isorp = 1, nssh(indna)
        call trescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, rhomx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            rho_off(imu,inu,mneigh,iatom) = rho_off(imu,inu,mneigh,iatom) + &
            &                               rhomx(imu,inu)*Qin(isorp,ialp)
            rho_off(inu,imu,jneigh,jatom) = rho_off(imu,inu,mneigh,iatom)
          end do ! imu
        end do ! inu
      end do ! isorp
    end do ! ineigh
  end do ! ialp
end subroutine average_ca_rho
