subroutine build_zw_on_na (in1, iatom, bcxcx, xc)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: bcxcx, arho_on, arhoi_on, rho_on, rhoi_on, numorb_max
  use M_fdata, only: nssh, Qneutral, lssh, nsh_max, num_orb
  implicit none
  integer, intent (in) :: in1
  integer, intent (in) :: iatom
  real, intent (out), dimension (numorb_max, numorb_max) :: bcxcx
  real, intent (out) :: xc
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
  real dexc
  real d2exc 
  real dmuxc
  real d2muxc
  real exc
  real muxc
  real dexci
  real d2exci
  real dmuxci
  real d2muxci
  real exci
  real muxci
  real q_mu
  real, dimension (nsh_max,nsh_max) :: arho
  real, dimension (nsh_max,nsh_max) :: arhoi
  real, dimension (numorb_max, numorb_max) :: denx
  real, dimension (numorb_max, numorb_max) :: deni
  xc = 0.0d0
  bcxcx = 0.0d0
  do inu = 1, nssh(in1)
   do imu = 1, nssh(in1)
    arho(imu,inu) = arho_on(imu,inu,iatom)
    arhoi(imu,inu) = arhoi_on(imu,inu,iatom)
   end do
  end do
  do inu = 1, num_orb(in1)
   do imu = 1, num_orb(in1)
    denx(imu,inu) = rho_on(imu,inu,iatom) 
    deni(imu,inu) = rhoi_on(imu,inu,iatom)
   end do
  end do
  n1 = 0
  do issh = 1, nssh(in1)
   l1 = lssh(issh,in1)
   n1 = n1 + l1 + 1
   call cepal (arhoi(issh,issh), exci, muxci, dexci, d2exci, dmuxci,d2muxci)
   call cepal (arho(issh,issh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
   do ind1 = -l1, l1
    imu = n1 + ind1
    bcxcx(imu,imu) = muxc
    bcxcx(imu,imu) = bcxcx(imu,imu) + dmuxc*(denx(imu,imu) - arho(issh,issh) )
    bcxcx(imu,imu) = bcxcx(imu,imu) - muxci - dmuxci*(deni(imu,imu) - arhoi(issh,issh) )
   end do
   
   q_mu = Qneutral(issh,in1) / (2*l1 + 1)
   do ind1 = -l1, l1
    imu = n1 + ind1
    xc = xc + q_mu*(exc - muxc + (dexc - dmuxc)*(denx(imu,imu) - arho(issh,issh)))
    xc = xc + q_mu*(muxci - exci - dexci*deni(imu,imu) + dexci*arhoi(issh,issh) + dmuxci*deni(imu,imu) - dmuxci*arhoi(issh,issh))
   end do
   n1 = n1 + l1
  end do 
  n1 = 0
  do issh = 1, nssh(in1)
   l1 = lssh(issh,in1)
   n1 = n1 + l1 + 1
   n2 = 0
   do jssh = 1, nssh(in1)
    l2 = lssh(jssh,in1)
    n2 = n2 + l2 + 1
    call cepal(arhoi(issh,jssh), exci, muxci, dexci, d2exci,dmuxci, d2muxci)
    call cepal(arho(issh,jssh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
    do ind1 = -l1, l1
     imu = n1 + ind1
     do ind2 = -l2, l2
      inu = n2 + ind2
      if (imu .ne. inu) then
       bcxcx(imu,inu) = dmuxc*denx(imu,inu) 
       bcxcx(imu,inu) = bcxcx(imu,inu) - dmuxci*deni(imu,inu) 
      end if
     end do 
    end do 
    n2 = n2 + l2
   end do 
   n1 = n1 + l1
  end do 
  return
end subroutine build_zw_on_na
