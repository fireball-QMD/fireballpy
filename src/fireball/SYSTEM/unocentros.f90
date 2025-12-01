!subroutine unocentros (in1, iatom, exc_1c, muexc_1c, dccexc_1c, mu1xc)
subroutine unocentros (in1, iatom, dccexc_1c, mu1xc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: degelec, getmssh, getlssh, getissh, numorb_max, Qin
  use M_fdata, only: nsh_max, nssh, num_orb, Qneutral, exc_1c_0, vxc_1c_0, gxc_1c, fxc_1c
!borrar nuxc1c, dnuxc1c, exc1c0, dexc1c, 
  implicit none
  integer, intent(in) :: iatom
  integer, intent(in) :: in1
!  real(double), intent(out) :: exc_1c  ! XC energy term of DCC
!  real(double), intent(out) :: muexc_1c      ! XC potential term of DCC
  real(double), intent(out) :: dccexc_1c     ! XC DCC term
  real(double), intent(out), dimension (numorb_max, numorb_max) :: mu1xc
  integer imu
  integer issh
  integer l1, l2
  integer m1, m2
  integer inu
  integer jssh
  integer kssh
!  real(double), dimension (nsh_max) :: dqi
  real(double), dimension (nsh_max) :: V_aux
  real(double), dimension (nsh_max) :: E_aux

!  exc_1c = 0.0d0
  dccexc_1c = 0.0d0
!  muexc_1c = 0.0d0
  mu1xc = 0.0d0
  !dqi = 0.0d0
  !do issh = 1, nssh(in1)
    !dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
    !dqi(issh) =0.0d0 !assemble_xczw
  !end do

  print*,'---- unocentros --',in1, iatom, num_orb(in1)

  do imu = 1,num_orb(in1)
    m1   = getmssh(degelec(iatom)+imu)
    l1   = getlssh(degelec(iatom)+imu)
    issh = getissh(degelec(iatom)+imu)
    do inu = 1,num_orb(in1)
      m2   = getmssh(degelec(iatom)+inu) 
      l2   = getlssh(degelec(iatom)+inu)
      jssh = getissh(degelec(iatom)+inu)
      if( m1 .eq. m2 .and. l1 .eq. l2 ) then
        !mu1xc(inu,imu) = nuxc1c(in1,jssh,issh)
        mu1xc(inu,imu) = vxc_1c_0(in1,jssh,issh)  
        do kssh = 1,nssh(in1)
          !mu1xc(inu,imu) = mu1xc(inu,imu) +  dnuxc1c(in1,jssh,issh,kssh)*dqi(kssh)
          mu1xc(inu,imu) = mu1xc(inu,imu) + Qin(kssh,iatom) * gxc_1c(in1,jssh,issh,kssh)
        enddo 
      endif
    end do
  end do
  v_aux = 0.d0
  do issh = 1,nssh(in1)
    v_aux(issh) = vxc_1c_0(in1,issh,issh)
    do kssh = 1,nssh(in1)
      v_aux(issh) = v_aux(issh) + Qin(kssh,iatom) * gxc_1c(in1,issh,issh,kssh)
    enddo
  enddo
  e_aux = 0.d0  
  do issh = 1,nssh(in1)
    e_aux(issh) = exc_1c_0(in1,issh) 
    do kssh = 1,nssh(in1)
      e_aux(issh) = e_aux(issh) + Qin(kssh,iatom) * fxc_1c(in1,issh,kssh)
    enddo
  enddo

  do issh = 1,nssh(in1)
    dccexc_1c = dccexc_1c + (e_aux(issh) - v_aux(issh))*Qin(issh,iatom)    
  enddo
   ! exc_1c = exc_1c + Qin(issh,iatom)*e_aux(issh)
   ! exc_1c = exc_1c + exc1c0(in1,issh,issh)*Qin(issh,iatom) 
   ! muexc_1c = muexc_1c + nuxc1c(in1,issh,issh)*Qin(issh,iatom)
    !dccexc_1c = dccexc_1c +  (exc1c0(in1,issh,issh) - nuxc1c(in1,issh,issh))*Qin(issh,iatom)
!    do jssh = 1,nssh(in1)
!      exc_1c = exc_1c + dexc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom) 
!      muexc_1c = muexc_1c +  dnuxc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom)
!      dccexc_1c = dccexc_1c +  ( dexc1c(in1,issh,issh,jssh) - dnuxc1c(in1,issh,issh,jssh) )*dqi(jssh)*Qin(issh,iatom)
!    end do
!  end do
  return
end subroutine unocentros
