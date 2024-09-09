subroutine unocentros (in1, iatom, exc_1c, muexc_1c, dccexc_1c, mu1xc)
  use iso_c_binding
  use M_system
  use M_fdata
  implicit none
  integer(c_long), intent(in) :: iatom
  integer(c_long), intent(in) :: in1
  real(c_double), intent(out) :: exc_1c  ! XC energy term of DCC
  real(c_double), intent(out) :: muexc_1c      ! XC potential term of DCC
  real(c_double), intent(out) :: dccexc_1c     ! XC DCC term
  real(c_double), intent(out), dimension (numorb_max, numorb_max) :: mu1xc
  integer(c_long) imu
  integer(c_long) issh
  integer(c_long) l1, l2
  integer(c_long) m1, m2
  integer(c_long) inu
  integer(c_long) jssh
  integer(c_long) kssh
  real(c_double), dimension (nsh_max) :: dqi

  exc_1c = 0.0d0
  dccexc_1c = 0.0d0
  muexc_1c = 0.0d0
  mu1xc = 0.0d0
  dqi = 0.0d0
  do issh = 1, nssh(in1)
    dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
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
        mu1xc(inu,imu) = nuxc1c(in1,jssh,issh)
        do kssh = 1,nssh(in1)
          mu1xc(inu,imu) = mu1xc(inu,imu) +  dnuxc1c(in1,jssh,issh,kssh)*dqi(kssh)
        enddo
      endif
    end do
  end do
  do issh = 1,nssh(in1)
    exc_1c = exc_1c + exc1c0(in1,issh,issh)*Qin(issh,iatom) 
    muexc_1c = muexc_1c + nuxc1c(in1,issh,issh)*Qin(issh,iatom)
    dccexc_1c = dccexc_1c +  (exc1c0(in1,issh,issh) - nuxc1c(in1,issh,issh))*Qin(issh,iatom)
    do jssh = 1,nssh(in1)
      exc_1c = exc_1c + dexc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom) 
      muexc_1c = muexc_1c +  dnuxc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom)
      dccexc_1c = dccexc_1c +  ( dexc1c(in1,issh,issh,jssh) - dnuxc1c(in1,issh,issh,jssh) )*dqi(jssh)*Qin(issh,iatom)
    end do
  end do
  return
end subroutine unocentros
