! This routine assembles all of the one-center exchange-interactions. The results are stored in vxc_1c and etotxc_1c
subroutine assemble_olsxc_1c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, etotxc_1c, Uexc_1c, Umuxc_1c, vxc_1c, neigh_self, numorb_max
  use M_fdata, only: num_orb
  implicit none
  integer iatom
  integer imu,inu
  integer in1
  integer matom
  real(double) dccexc_1c
  real(double) exc_1c
  real(double) muexc_1c
  real(double), dimension (numorb_max,numorb_max) :: mu1xc
  etotxc_1c = 0.0d0
  Uexc_1c = 0.0d0
  Umuxc_1c = 0.0d0
  vxc_1c = 0.0d0
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    in1 = imass(iatom)
    call unocentros (in1, iatom, exc_1c, muexc_1c,  dccexc_1c, mu1xc)
    etotxc_1c = etotxc_1c + dccexc_1c
    do imu = 1, num_orb(in1)
      do inu = 1, num_orb(in1)
        vxc_1c(imu,inu,matom,iatom) =  vxc_1c(imu,inu,matom,iatom) + mu1xc(imu,inu) 
      end do
    end do
  end do
end subroutine assemble_olsxc_1c
