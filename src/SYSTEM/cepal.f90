! This routine calculates the lda exchange-correlation energy, potential and the derivative of the potential.  The form of the
! functional is that of Ceperley-Alder as parameterized by Perdew-Zunger.
! Phys. Rev. B23, 5048 (1981). Units are Hartree a.u. (but see below)
subroutine cepal (rh, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
  use M_constants
  implicit none
  real(8), intent (in) :: rh
  real(8), intent (out) :: dexc
  real(8), intent (out) :: d2exc
  real(8), intent (out) :: dmuxc
  real(8), intent (out) :: exc
  real(8), intent (out) :: muxc
  real(8), intent (out) :: d2muxc
  real(8), parameter :: eps = 1.0d-3
  real(8), parameter :: delta_rh = 1.0d-6
  real(8) d2nec
  real(8) d2nex
  real(8) d3nec
  real(8) d3nex
  real(8) dec
  real(8) ddec
  real(8) d2dec
  real(8) den
  real(8) dden
  real(8) d2den
  real(8) d3den
  real(8) ex
  real(8) rho_third
  real(8) rho
  real(8) rhx
  real(8) rs
  real(8) rsl
  real(8) sqrs
  real(8) hartree1
 
  exc = 0.0d0
  muxc = 0.0d0
  dexc = 0.0d0
  dmuxc = 0.0d0
  d2muxc = 0.0d0
  d2exc = 0.0d0

  rhx = sqrt(rh*rh + delta_rh)
  rho = rhx*(abohr**3)
  rho_third=rho**(1.0e0/3.0e0)
  rs = 0.62035049d0/rho_third
  if (rho .lt. 0.23873241d0) then
   sqrs = sqrt(rs)
   den = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs    ! Effective density
   exc = -0.4581652d0/rs - 0.1423d0/den
   muxc = exc - rs*(0.15273333d0/rs**2 + (0.02497128d0/sqrs + 0.01581427d0)/den**2)
   dden = 1.0529d0/(2.0d0*sqrs) + 0.3334d0
   d2den = (-0.5d0)*1.0529d0/(2.0d0*rs*sqrs)
   d3den = (0.75d0)*1.0529d0/(2.0d0*rs*rs*sqrs)
   dec = 0.1423d0*dden/(den*den)
   ddec = -2.0d0*0.1423d0*dden*dden/(den**3) + 0.1423d0*d2den/(den*den)
   d2dec = 6.0d0*0.1423d0*(dden*3)/(den**4) - 6.0d0*0.1423d0*dden*d2den/(den**3) + 0.1423d0*d3den/(den*den) 
  else
   rsl = log(rs)
   exc = -0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs + 0.002d0*rs*rsl
   muxc = exc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))
   dec = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
   ddec = -0.0311d0/(rs*rs) + 0.0020d0/rs
   d2dec = 2.0d0*0.0311d0/(rs*rs*rs) - 0.0020d0/(rs*rs)
  end if
  ex = -0.7385587664d0*rho_third 
  dexc = (muxc - exc)/rho
  d2nec = (4.0d0*rs/(9.0d0*rho*rho))*dec + (rs*rs/(9.0d0*rho*rho))*ddec
  d2nex = -(2.0d0/(9.0d0*rho*rho))*ex
  dmuxc = 2.0d0*dexc + rho*(d2nex + d2nec)
  d3nec = (-28.0d0*rs/(27.0d0*rho*rho*rho))*dec + (-4.0d0*rs*rs/(9.0d0*rho*rho*rho))*ddec + (rs*rs*rs/(-27.0d0*rho*rho*rho))*d2dec
  d3nex = (10.0d0/(27.0d0*rho*rho*rho))*ex
  d2muxc = 3.0*(d2nex + d2nec) + rho*(d3nex + d3nec)
  d2exc = d2nex + d2nec
  hartree1 = eq2/abohr
  exc = exc*hartree1
  muxc = muxc*hartree1
  dexc = dexc*hartree1*(abohr)**3
  d2exc = d2exc*hartree1*(abohr)**6
  dmuxc = dmuxc*hartree1*(abohr)**3
  d2muxc = d2muxc*hartree1*(abohr)**6
  return
end subroutine cepal
