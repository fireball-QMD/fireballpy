! This routine assembles all of the one-center exchange-interactions. The results are stored in vxc_1c and etotxc_1c
subroutine assemble_olsxc_1c ()
  use M_system
  use M_fdata, only: num_orb
  implicit none
  integer iatom
  integer imu,inu
  integer in1
  integer matom
  real(8) dccexc_1c
  real(8) exc_1c
  real(8) muexc_1c
  real(8), dimension (numorb_max,numorb_max) :: mu1xc
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

