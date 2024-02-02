subroutine buildh (itestrange,testrange)
  use M_system
  implicit none
  integer, intent (in) :: itestrange
  real, intent (in) :: testrange
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
        end do ! do imu
      end do ! do inu
    end do ! do ineigh
    if (V_intra_dip .eq. 1) then
      matom = neigh_self(iatom)
      do inu = 1, num_orb(in1)
        do imu = 1, num_orb(in1)
          h_mat(imu,inu,matom,iatom) = h_mat(imu,inu,matom,iatom) + Vdip_1c(imu,inu,iatom)
        end do ! do imu
      end do ! do inu
    end if ! if (V_intra_dip .eq. 1) 
  end do ! do iatom
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      if (itestrange .eq. 0) then
        distance = sqrt((ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2 + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2  + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2)
        if (distance .gt. testrange) then
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              h_mat(imu,inu,ineigh,iatom) = 0.0d0
              t_mat(imu,inu,ineigh,iatom) = 0.0d0
              s_mat(imu,inu,ineigh,iatom) = 0.0d0
              vna(imu,inu,ineigh,iatom) = 0.0d0
              vxc(imu,inu,ineigh,iatom) = 0.0d0
              ewaldlr(imu,inu,ineigh,iatom) = 0.0d0
              ewaldsr(imu,inu,ineigh,iatom) = 0.0d0
              ewaldqmmm(imu,inu,ineigh,iatom) = 0.0d0
            end do
          end do
        end if
      end if
    end do ! do ineigh
  end do ! do iatom
  return
end subroutine buildh

