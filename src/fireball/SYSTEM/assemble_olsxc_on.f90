subroutine assemble_olsxc_on ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: natoms, imass, uxcdcc_ols, neigh_self, numorb_max, vxc
  use M_fdata, only: num_orb
  implicit none
  integer iatom, matom
  integer imu
  integer in1, in3
  integer inu
  real(dp), dimension (numorb_max, numorb_max) :: bcxcx
  real(dp) xc
  vxc = 0.0_dp
  uxcdcc_ols = 0.0_dp
  bcxcx  = 0.0_dp
  xc = 0.0_dp
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    in1 = imass(iatom)
    call build_ca_olsxc_on (in1, iatom, bcxcx, xc)
    uxcdcc_ols = uxcdcc_ols + xc
    in3 = in1
    do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
        vxc(imu,inu,matom,iatom) =  vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
      end do
    end do
  end do 
end subroutine assemble_olsxc_on
