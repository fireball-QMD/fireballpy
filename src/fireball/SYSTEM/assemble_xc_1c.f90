subroutine assemble_xc_1c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: natoms, neigh_self, imass, numorb_max, &
                    etotxc_1c, Uexc_1c, Umuxc_1c, vxc_1c,g_xc
  use M_fdata, only: num_orb
  implicit none
 
  integer iatom
  integer imu,inu
  integer in1
  integer matom
  integer ixc
  real(double) dccexc_1c
  real(double) exc_1c
  real(double) muexc_1c
  real(double), dimension (numorb_max,numorb_max) :: mu1xc
  etotxc_1c = 0.0d0
  Uexc_1c = 0.0d0
  Umuxc_1c = 0.0d0
  vxc_1c = 0.0d0
  g_xc = 0.0d0
  do iatom = 1, natoms
   matom = neigh_self(iatom)
   in1 = imass(iatom)
   ixc = 4
   print*,'assemble_xc_1c iatom =', iatom, in1
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
 
