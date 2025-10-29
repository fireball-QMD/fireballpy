subroutine assemble_zw_on_na ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, natoms, neigh_self, imass, vxc, uxcdcc_zw
  use M_fdata, only: num_orb
  implicit none
  integer iatom
  integer imu
  integer in1, in3
  integer inu
  integer matom
  real, dimension (numorb_max, numorb_max) :: bcxcx
  real xc
  vxc = 0.0d0
  uxcdcc_zw = 0.0d0   !this quantity is initialized here
  bcxcx  = 0.0d0
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    in1 = imass(iatom)
    call build_zw_on_na (in1, iatom, bcxcx, xc)
    ! double-counting xc correction
    uxcdcc_zw = uxcdcc_zw + xc
    in3 = in1
    do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
        vxc(imu,inu,matom,iatom) =  vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
      end do
    end do
  end do ! End loop over iatom.
  return
end subroutine assemble_zw_on_na
