subroutine buildh ()
  use iso_c_binding
  use M_system
  use m_fdata, only: num_orb 
  implicit none
  integer(c_long) katom
  integer(c_long) iatom
  integer(c_long) ierror
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) ineigh
  integer(c_long) matom
  integer(c_long) inu
  integer(c_long) jatom
  integer(c_long) mbeta
  real(c_double) distance
  real(c_double), dimension (numorb_max, numorb_max) :: htemp
  real(c_double), dimension (numorb_max, numorb_max) :: stemp
  real(c_double), dimension (3) :: dvec
  integer(c_long) issh
  integer(c_long) numorb
  integer(c_long) jatom0
  integer(c_long) ineigh0
  integer(c_long) mbeta0
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
          h_mat(imu,inu,ineigh,iatom) = h_mat(imu,inu,ineigh,iatom)  + vca(imu,inu,ineigh,iatom) + vxc_ca(imu,inu,ineigh,iatom) + ewaldlr(imu,inu,ineigh,iatom) - ewaldsr(imu,inu,ineigh,iatom) + ewaldqmmm(imu,inu,ineigh,iatom)
          !AQUI  ewaldqmmm
        end do 
      end do 
    end do 
  end do

 
end subroutine buildh

