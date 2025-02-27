subroutine assemble_qmmm ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, eqmmm, dip, neigh_b, neigh_j, neighn, s_mat, ewaldqmmm, xl, qmmm_qm_natoms, &
    & qmmm_qm_xcrd
  use M_fdata, only: nssh, Qneutral, num_orb
  implicit none
  integer(c_long) :: iatom, imu, inu, in1, in2, in3, ineigh, issh, jatom, katom, mbeta
  real(c_double) :: distance12, dij, dterm, sterm, dq3, v
  real(c_double), dimension (3) :: r1, r2
  real(c_double), dimension (natoms) :: sub_ewaldqmmm
  real(c_double), external :: sf

  ewaldqmmm = 0.0d0
  sub_ewaldqmmm = 0.0d0
  do iatom = 1, natoms
    do katom = 1, qmmm_qm_natoms
      dij = norm2(qmmm_qm_xcrd(:,katom) - ratom(:,iatom))
      v = qmmm_qm_xcrd(4, katom) * sf(dij)
      sub_ewaldqmmm(iatom) = sub_ewaldqmmm(iatom) - v
    end do
  end do
  do iatom = 1,natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      mbeta = neigh_b(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      distance12 = norm2(r2 - r1)
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
          if (distance12 .gt. 1.0d-4) then
            dterm = dip(imu,inu,ineigh,iatom)/distance12
          else
            dterm = 0.0d0
          end if
          ewaldqmmm(imu,inu,ineigh,iatom) =  ewaldqmmm(imu,inu,ineigh,iatom)  + (sterm - dterm)*sub_ewaldqmmm(iatom)*eq2  + (sterm + dterm)*sub_ewaldqmmm(jatom)*eq2
        end do
      end do
    end do
  end do
  eqmmm = 0.0d0
  do iatom = 1,natoms
    in3 = imass(iatom)
    dq3 = 0.0d0
    do issh = 1, nssh(in3)
      dq3 = dq3  + Qneutral(issh,in3)
    end do
    eqmmm = eqmmm - dq3*sub_ewaldqmmm(iatom)*eq2
  end do
end subroutine assemble_qmmm
