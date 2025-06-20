subroutine geth (overlap, hamiltonian)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neighn, neigh_j, h_mat, s_mat, neighPP_j, neighn, neighPPn, &
    & vnl, degelec, norbitals
  use M_fdata, only: num_orb
  implicit none
  integer :: iatom, jatom, imu, inu, in1, in2, ineigh, jmu, jnu
  real(double), dimension (norbitals, norbitals), intent(out) :: overlap, hamiltonian

  overlap = 0.0d0
  hamiltonian = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          overlap(jmu,jnu) = overlap(jmu,jnu) + s_mat(imu,inu,ineigh,iatom)
          hamiltonian(jmu,jnu) = hamiltonian(jmu,jnu) + h_mat(imu,inu,ineigh,iatom)
        end do
      end do
    end do

    do ineigh = 1, neighPPn(iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          hamiltonian(jmu,jnu) = hamiltonian(jmu,jnu) + vnl(imu,inu,ineigh,iatom)
        end do
      end do
    end do
  end do
end subroutine geth
