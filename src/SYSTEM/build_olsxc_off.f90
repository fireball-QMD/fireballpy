subroutine build_olsxc_off (in1, in2, den1x, denx, sx, ineigh,iatom, bcxcx)
  use M_constants, only: wp
  use M_system
  use M_fdata, only: nssh,lssh,nsh_max
  implicit none
  integer, intent (in) :: in1
  integer, intent (in) :: in2
  integer, intent (in) :: ineigh
  integer, intent (in) :: iatom
  real(wp), intent (in), dimension (numorb_max, numorb_max) :: denx
  real(wp), intent (in), dimension (numorb_max, numorb_max) :: den1x
  real(wp), intent (in), dimension (numorb_max, numorb_max) :: sx
  real(wp), intent (out), dimension (numorb_max, numorb_max) :: bcxcx
  integer imu
  integer ind1
  integer ind2
  integer inu
  integer issh
  integer jssh
  integer l1
  integer l2
  integer n1
  integer n2
  real(wp) dexc
  real(wp) d2exc
  real(wp) dmuxc
  real(wp) d2muxc
  real(wp) exc
  real(wp) muxc
  real(wp) dmuxcij
  real(wp) muxcij
  real(wp), dimension (nsh_max, nsh_max) :: dens
  real(wp), dimension (nsh_max, nsh_max) :: densij
  bcxcx = 0.0d0
  do issh = 1, nssh(in1)
    do jssh = 1, nssh(in2)
      dens(issh,jssh) = arho_off(issh,jssh,ineigh,iatom)
      densij(issh,jssh) = arhoij_off(issh,jssh,ineigh,iatom)
    enddo
  enddo
  n1 = 0
  do issh = 1, nssh(in1)
    l1 = lssh(issh,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do jssh = 1, nssh(in2)
      l2 = lssh(jssh,in2)
      n2 = n2 + l2 + 1
      call cepal (dens(issh,jssh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
      call cepal (densij(issh,jssh), exc, muxcij, dexc, d2exc, dmuxcij, d2muxc)
      do ind1 = -l1, l1
        imu = n1 + ind1
        do ind2 = -l2, l2
          inu = n2 + ind2
          bcxcx(imu,inu) = muxc*sx(imu,inu)
          bcxcx(imu,inu) = bcxcx(imu,inu) + dmuxc*(denx(imu,inu) - dens(issh,jssh)*sx(imu,inu))
          bcxcx(imu,inu) = bcxcx(imu,inu) - muxcij*sx(imu,inu) - dmuxcij*(den1x(imu,inu) - densij(issh,jssh)*sx(imu,inu))
        end do ! do ind2 = -l2, l2
      end do ! do ind1 = -l1, l1
      n2 = n2 + l2
    end do ! do jssh = 1, nssh(in2)
    n1 = n1 + l1
  end do !do issh = 1, nssh(in1)
  return
end subroutine build_olsxc_off

