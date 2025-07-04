subroutine unocentros (in1, iatom, exc_1c, muexc_1c, dccexc_1c, mu1xc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: degelec, getmssh, getlssh, getissh, numorb_max, Qin
  use M_fdata, only: nsh_max, nssh, num_orb, nuxc1c, dnuxc1c, exc1c0, dexc1c, Qneutral
  implicit none
  integer, intent(in) :: iatom
  integer, intent(in) :: in1
  real(double), intent(out) :: exc_1c  ! XC energy term of DCC
  real(double), intent(out) :: muexc_1c      ! XC potential term of DCC
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
