subroutine buildh ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, vxc_1c, ewaldsr, neigh_j, neighn, vxc, vxc_ca, vca, ewaldlr, h_mat, t_mat, vna, ewaldqmmm
  use M_fdata, only: num_orb
  implicit none
  integer :: iatom, jatom, imu, inu, ineigh, in1, in2

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
            & vxc_ca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom) + &
            & ewaldqmmm(imu,inu,ineigh,iatom)
        end do
      end do
    end do
  end do
end subroutine buildh
