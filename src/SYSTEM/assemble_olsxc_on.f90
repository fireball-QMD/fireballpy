! ===========================================================================
! This routine assembles the two-center exchange-correlation
! (on-site - atom case) for the average density approximation. 
! This subroutine could be easily incorporated in assemble_2c.f90 (JOM)
! The double-counting xc (uxcdcc) is also calculated here.
subroutine assemble_olsxc_on (uxcdcc)
  use M_system
  implicit none
  real, intent (out) :: uxcdcc
  integer iatom
  integer imu
  integer in1, in3
  integer inu
  real, dimension (numorb_max, numorb_max) :: bcxcx
  real xc
  vxc = 0.0d0
  if (itheory .eq. 1) vxc_ca = 0.0d0
  uxcdcc = 0.0d0
  bcxcx  = 0.0d0
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    in1 = imass(iatom)
    call build_ca_olsxc_on (in1, iatom, bcxcx, xc)
    uxcdcc = uxcdcc + xc
    in3 = in1
    do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
        vxc_ca(imu,inu,matom,iatom) =  vxc_ca(imu,inu,matom,iatom) + bcxcx(imu,inu)
      end do
    end do
  end do 
  return
end subroutine assemble_olsxc_on
