subroutine Dassemble_qmmm ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, dip, dipp, neigh_b, neigh_j, neighn, sp_mat, rho, s_mat, &
    & xl, flrew_qmmm, qmmm_qm_natoms, qmmm_dxyzcl, qmmm_qm_xcrd
  use M_fdata, only  : nssh, Qneutral, num_orb
  implicit none
  integer(c_long) :: iatom, imu, in1, in2, in3, ineigh, inu, issh, jatom, jmu, katom, mbeta
  real(c_double) :: distance12, dij, dq3, dterm, sterm, sff
  real(c_double), dimension (3) :: dewaldlr_i_qmmm, dewaldlr_j_qmmm, dpterm, r1, r2, rhat12, spterm, vij, u
  real(c_double), dimension (natoms) :: sub_ewaldqmmm
  real(c_double), dimension (3, natoms) :: sub_dewaldqmmm
  real(c_double), dimension (3, natoms, qmmm_qm_natoms) :: dqmmm
  real(c_double), external :: sf

  flrew_qmmm = 0.0d0
  sub_ewaldqmmm = 0.0d0
  sub_dewaldqmmm = 0.0d0
  dqmmm(:,:,:) = 0.0d0
  qmmm_dxyzcl = 0.0d0
  do iatom = 1, natoms
    do katom = 1, qmmm_qm_natoms
      vij = qmmm_qm_xcrd(:, katom) - ratom(:, iatom)
      dij = norm2(vij)
      u = vij/dij
      sff = sf(dij)
      sub_ewaldqmmm(iatom) = sub_ewaldqmmm(iatom) - qmmm_qm_xcrd(4,katom)*sff
      do jmu = 1, 3
        dqmmm(jmu,iatom,katom) = -(u(jmu)*sff*sff)
        sub_dewaldqmmm(jmu,iatom) = sub_dewaldqmmm(jmu,iatom) - qmmm_qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)
      end do
    end do
  end do
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      mbeta = neigh_b(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      distance12 = norm2(r2 - r1)
      if (distance12 .gt. 1.0d-4) rhat12(:) = (r2(:) - r1(:))/distance12
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
          if (distance12 .gt. 1.0d-4) then
            dterm = dip(imu,inu,ineigh,iatom)/distance12
            dpterm(:) = dipp(:,imu,inu,ineigh,iatom)/distance12 + dip(imu,inu,ineigh,iatom)*rhat12(:)/distance12**2
            spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)/2.0d0
            dewaldlr_i_qmmm(:) = (sterm - dterm)*sub_dewaldqmmm(:,iatom) &
              & + (spterm(:) - dpterm(:))*sub_ewaldqmmm(iatom) &
              & + (spterm(:) + dpterm(:))*sub_ewaldqmmm(jatom)
            dewaldlr_j_qmmm(:) = - (spterm(:) - dpterm(:))*sub_ewaldqmmm(iatom) &
              & + (sterm + dterm)*sub_dewaldqmmm(:,jatom) &
              & - (spterm(:) + dpterm(:))*sub_ewaldqmmm(jatom)
          else
            dterm = 0.0d0
            dpterm(:) = 0.0d0
            spterm(:) = 0.0d0
            dewaldlr_i_qmmm(:) = sterm*sub_dewaldqmmm(:,iatom)
            dewaldlr_j_qmmm(:) = sterm*sub_dewaldqmmm(:,jatom)
          end if
          do jmu = 1, 3
            flrew_qmmm(jmu,iatom) = flrew_qmmm(jmu,iatom) - rho(imu,inu,ineigh,iatom)*dewaldlr_i_qmmm(jmu)*eq2
            flrew_qmmm(jmu,jatom) = flrew_qmmm(jmu,jatom) - rho(imu,inu,ineigh,iatom)*dewaldlr_j_qmmm(jmu)*eq2
          end do
          do katom = 1, qmmm_qm_natoms
            do jmu = 1, 3
              qmmm_dxyzcl(jmu,katom) = qmmm_dxyzcl(jmu,katom) &
                & - ((sterm-dterm)*dqmmm(jmu,iatom,katom) &
                & +(sterm+dterm)*dqmmm(jmu,jatom,katom))*rho(imu,inu,ineigh,iatom)*qmmm_qm_xcrd(4,katom)*eq2
            end do
          end do
        end do
      end do
    end do
  end do
  do iatom = 1, natoms
    in3 = imass(iatom)
    dq3 = 0.0d0
    do issh = 1, nssh(in3)
      dq3 = dq3  + Qneutral(issh,in3)
    end do
    do katom = 1, qmmm_qm_natoms
      do jmu = 1, 3
        qmmm_dxyzcl(jmu,katom) = qmmm_dxyzcl(jmu,katom) + dq3*qmmm_qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2
        flrew_qmmm(jmu,iatom) = flrew_qmmm(jmu,iatom) + dq3*qmmm_qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2
      end do
    end do
  end do
end subroutine Dassemble_qmmm
