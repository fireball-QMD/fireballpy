subroutine buildh ()
  use M_system
  use m_fdata, only: num_orb, V_intra_dip 
  implicit none
  integer katom
  integer iatom
  integer iatomstart
  integer ierror
  integer imu
  integer in1
  integer in2
  integer ineigh
  integer matom
  integer inu
  integer jatom
  integer mbeta
  integer my_proc
  integer natomsp
  real distance
  real, dimension (numorb_max, numorb_max) :: htemp
  real, dimension (numorb_max, numorb_max) :: stemp
  real, dimension (3) :: dvec
  integer issh
  integer numorb
  integer jatom0
  integer ineigh0
  integer mbeta0
  write(*,*)'XXX in buildh vna',vna
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          h_mat(imu,inu,ineigh,iatom) =  t_mat(imu,inu,ineigh,iatom) + vna(imu,inu,ineigh,iatom)  + vxc(imu,inu,ineigh,iatom) + vxc_1c(imu,inu,ineigh,iatom)
        end do
      end do
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          !AQUI  ewaldqmmm
          h_mat(imu,inu,ineigh,iatom) = h_mat(imu,inu,ineigh,iatom)  + vca(imu,inu,ineigh,iatom) + vxc_ca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom) + ewaldqmmm(imu,inu,ineigh,iatom)
        end do ! do imu
      end do ! do inu
    end do ! do ineigh
  end do ! do iatom
  write(*,*)'XXX buildh vna',vna
  write(*,*)'XXX buildh vxc',vxc
  write(*,*)'XXX buildh vxc_ca',vxc_ca
  write(*,*)'XXX buildh vxc_1c',vxc_1c
  write(*,*)'XXX buildh t_mat',t_mat
  write(*,*)'XXX buildh h_mat',h_mat
end subroutine buildh

