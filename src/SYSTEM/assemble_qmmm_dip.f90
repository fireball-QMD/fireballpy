subroutine assemble_qmmm_dip () 
  use M_constants, only : eq2, wp
  use M_system
  use M_fdata, only  : nssh, Qneutral, num_orb
  implicit none
  integer :: iatom, imu, inu, in1, in2, in3, ineigh, issh, jatom, katom, mbeta
  real(wp) :: dij, dterm, sterm, dq3, dq4, x, emnpl
  real(wp), dimension (3) :: rna, rnabc, r1, r2, r21

  ewaldqmmm = 0.0d0
  emnpl = 0.0d0
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) !+ xl(:,mbeta)
      in2 = imass(jatom)
      do katom = 1, qmmm_qm_mm_pairs
        rna(1) = qmmm_qm_xcrd(1,katom)
        rna(2) = qmmm_qm_xcrd(2,katom)
        rna(3) = qmmm_qm_xcrd(3,katom)
        dq3 = - qmmm_qm_xcrd(4,katom) ! charge in amber have opposite sign
        r21(:) = r2(:) - r1(:)
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            sterm = s_mat(imu,inu,ineigh,iatom)
            dterm = dipc(1,imu,inu,ineigh,iatom)*rnabc(1) + dipc(2,imu,inu,ineigh,iatom)*rnabc(2) + dipc(3,imu,inu,ineigh,iatom)*rnabc(3)
            emnpl = dq3*sterm/x + dq3*dterm/(x*x*x)
            ewaldqmmm(imu,inu,ineigh,iatom) = ewaldqmmm(imu,inu,ineigh,iatom) + emnpl*eq2
          end do !end do imu = 1, num_orb(in1)
        end do  ! end do inu = 1, num_orb(in2)
      end do    ! end do katom: mm atom
    end do  ! end do ineigh = 1, neighn(iatom)
  end do   ! end do iatom = 1, natoms
  eqmmm = 0.0d0
  do iatom = 1, natoms
    in3 = imass(iatom)
    dq4 = 0.0d0
    do issh = 1, nssh(in3)
      dq4 = dq4  + Qneutral(issh,in3)
    end do
    do katom = 1, qmmm_qm_mm_pairs
      dij = sqrt( (qmmm_qm_xcrd(1,katom)-ratom(1,iatom))*(qmmm_qm_xcrd(1,katom)-ratom(1,iatom)) &
        & + (qmmm_qm_xcrd(2,katom)-ratom(2,iatom))*(qmmm_qm_xcrd(2,katom)-ratom(2,iatom)) &
        & + (qmmm_qm_xcrd(3,katom)-ratom(3,iatom))*(qmmm_qm_xcrd(3,katom)-ratom(3,iatom)) )
      eqmmm = eqmmm + (dq4*qmmm_qm_xcrd(4,katom) / dij)*eq2
    end do
  end do
end subroutine assemble_qmmm_dip
