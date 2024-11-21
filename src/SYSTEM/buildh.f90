subroutine buildh ()
  use iso_c_binding
  use M_system, only: natoms, imass, vxc_1c, ewaldsr, neigh_j, neighn, vxc, vxc_ca, vca, ewaldlr, h_mat, t_mat, vna, ewaldqmmm
  use m_fdata, only: num_orb
  implicit none
  integer(c_long) iatom
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) ineigh
  integer(c_long) inu
  integer(c_long) jatom

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
