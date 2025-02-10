subroutine assemble_olsxc_on ()
  use iso_c_binding
  use M_system, only: natoms, imass, uxcdcc_ols, neigh_self, numorb_max, vxc, vxc_ca
  use M_fdata, only: num_orb
  implicit none
  integer(c_long) iatom, matom
  integer(c_long) imu
  integer(c_long) in1, in3
  integer(c_long) inu
  real(c_double), dimension (numorb_max, numorb_max) :: bcxcx
  real(c_double) xc
  vxc = 0.0d0
  vxc_ca = 0.0d0
  uxcdcc_ols = 0.0d0
  bcxcx  = 0.0d0
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    in1 = imass(iatom)
    call build_ca_olsxc_on (in1, iatom, bcxcx, xc)
    uxcdcc_ols = uxcdcc_ols + xc
    in3 = in1
    do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
        vxc_ca(imu,inu,matom,iatom) =  vxc_ca(imu,inu,matom,iatom) + bcxcx(imu,inu)
      end do
    end do
  end do 
end subroutine assemble_olsxc_on
