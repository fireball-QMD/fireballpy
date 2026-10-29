subroutine assemble_xc_1c ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: eq2
  use M_system, only: natoms, neigh_self, imass, numorb_max, &
                    etotxc_1c, Uexc_1c, Umuxc_1c, vxc_1c, g_xc, Kscf
  use M_fdata, only: num_orb
  implicit none
 
  integer iatom
  integer imu,inu
  integer in1
  integer matom
  real(dp) dccexc_1c
  real(dp) exc_1c
  real(dp) muexc_1c
  real(dp), dimension (numorb_max,numorb_max) :: mu1xc
  etotxc_1c = 0.0_dp
  Uexc_1c = 0.0_dp
  Umuxc_1c = 0.0_dp
  vxc_1c = 0.0_dp
  if (Kscf .eq. 1) g_xc = 0.0_dp
  do iatom = 1, natoms
   matom = neigh_self(iatom)
   in1 = imass(iatom)
   call unocentros (in1, iatom, dccexc_1c, mu1xc)
   etotxc_1c = etotxc_1c + dccexc_1c
   do imu = 1, num_orb(in1)
    do inu = 1, num_orb(in1)
     vxc_1c(imu,inu,matom,iatom) = vxc_1c(imu,inu,matom,iatom) + mu1xc(imu,inu) 
    end do
   end do
  end do
  return
end subroutine assemble_xc_1c
 
