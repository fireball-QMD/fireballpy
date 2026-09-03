!subroutine unocentros (in1, iatom, exc_1c, muexc_1c, dccexc_1c, mu1xc)
subroutine unocentros (in1, iatom, dccexc_1c, mu1xc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: degelec, getmssh, getlssh, getissh, numorb_max, Qin, neigh_self, g_xc, get_shell_ofatom_issh, Kscf
  use M_fdata, only: nsh_max, nssh, num_orb, Qref, exc_1c_0, vxc_1c_0, gxc_1c, fxc_1c
!borrar nuxc1c, dnuxc1c, exc1c0, dexc1c, 
  implicit none
  integer, intent(in) :: iatom
  integer, intent(in) :: in1
  real(double), intent(out) :: dccexc_1c     ! XC DCC term
  real(double), intent(out), dimension (numorb_max, numorb_max) :: mu1xc
  integer imu
  integer issh
  integer l1, l2
  integer m1, m2
  integer inu
  integer jssh
  integer kssh
  real(double), dimension (nsh_max) :: dqi

  dccexc_1c = 0.0d0
  mu1xc = 0.0d0
  dqi = 0.0d0
  do issh = 1, nssh(in1)
    dqi(issh) = (Qin(issh,iatom) - Qref(issh,in1))
  end do

  do imu = 1,num_orb(in1)
    m1   = getmssh(degelec(iatom)+imu)
    l1   = getlssh(degelec(iatom)+imu)
    issh = getissh(degelec(iatom)+imu)
    do inu = 1,num_orb(in1)
      m2   = getmssh(degelec(iatom)+inu) 
      l2   = getlssh(degelec(iatom)+inu)
      jssh = getissh(degelec(iatom)+inu)
      if( m1 .eq. m2 .and. l1 .eq. l2 ) then
         mu1xc(inu,imu) = vxc_1c_0(jssh,issh,in1)  
        do kssh = 1,nssh(in1)
          mu1xc(inu,imu) = mu1xc(inu,imu) + dqi(kssh) * gxc_1c(jssh,issh,kssh,in1)
          if (Kscf .eq. 1) g_xc(get_shell_ofatom_issh(iatom,kssh),inu,imu,neigh_self(iatom),iatom) = gxc_1c(jssh,issh,kssh,in1)
        enddo 
      endif
    end do
  end do
  do issh = 1, nssh(in1)
    dccexc_1c = dccexc_1c + (exc_1c_0(issh,in1) - vxc_1c_0(issh,issh,in1))*Qin(issh,iatom)
    do jssh = 1, nssh(in1)
      dccexc_1c = dccexc_1c + (fxc_1c(issh,jssh,in1) - gxc_1c(issh,issh,jssh,in1))*dqi(jssh)*Qin(issh,iatom)
    end do
  end do
end subroutine unocentros
