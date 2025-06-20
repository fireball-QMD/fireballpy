! This routine calculates the lda exchange-correlation energy, potential and the derivative of the potential.  The form of the
! functional is that of Ceperley-Alder as parameterized by Perdew-Zunger.
! Phys. Rev. B23, 5048 (1981). Units are Hartree a.u. (but see below)
subroutine cepal (rh, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2, abohr
  implicit none
  real(double), intent (in) :: rh
  real(double), intent (out) :: dexc
  real(double), intent (out) :: d2exc
  real(double), intent (out) :: dmuxc
  real(double), intent (out) :: exc
  real(double), intent (out) :: muxc
  real(double), intent (out) :: d2muxc
  real(double), parameter :: eps = 1.0d-3
  real(double), parameter :: delta_rh = 1.0d-6
  real(double) d2nec
  real(double) d2nex
  real(double) d3nec
  real(double) d3nex
  real(double) dec
  real(double) ddec
  real(double) d2dec
  real(double) den
  real(double) dden
  real(double) d2den
  real(double) d3den
  real(double) ex
  real(double) rho_third
  real(double) rho
  real(double) rhx
  real(double) rs
  real(double) rsl
  real(double) sqrs
  real(double) hartree1
 
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
