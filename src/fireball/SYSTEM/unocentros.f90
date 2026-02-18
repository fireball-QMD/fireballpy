!subroutine unocentros (in1, iatom, exc_1c, muexc_1c, dccexc_1c, mu1xc)
subroutine unocentros (in1, iatom, dccexc_1c, mu1xc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: degelec, getmssh, getlssh, getissh, numorb_max, Qin, neigh_self, g_xc
  use M_fdata, only: nsh_max, nssh, num_orb, Qneutral, exc_1c_0, vxc_1c_0, gxc_1c, fxc_1c
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
  real(double), dimension (nsh_max) :: V_aux
  real(double), dimension (nsh_max) :: E_aux

  dccexc_1c = 0.0d0
  mu1xc = 0.0d0
  dqi = 0.0d0

  do issh = 1, nssh(in1)
    dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
  end do

  print*,'---- unocentros --',in1, iatom

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
          g_xc(inu,imu,kssh,iatom,neigh_self(iatom),iatom) = &
          & g_xc(inu,imu,kssh,iatom,neigh_self(iatom),iatom) + gxc_1c(jssh,issh,kssh,in1)
        enddo 
      endif
    end do
  end do
  V_aux = 0.d0
  do issh = 1,nssh(in1)
    V_aux(issh) = vxc_1c_0(issh,issh,in1)
    do kssh = 1,nssh(in1)
      V_aux(issh) = V_aux(issh) + dqi(kssh) * gxc_1c(issh,issh,kssh,in1)
    enddo
  enddo
  e_aux = 0.d0  
  do issh = 1,nssh(in1)
    E_aux(issh) = exc_1c_0(issh,in1) 
    do kssh = 1,nssh(in1)
      E_aux(issh) = E_aux(issh) + dqi(kssh) * fxc_1c(issh,kssh,in1)
    enddo
  enddo

  do issh = 1,nssh(in1)
    dccexc_1c = dccexc_1c + (E_aux(issh) - V_aux(issh))*Qin(issh,iatom)    
  enddo
  return
end subroutine unocentros
