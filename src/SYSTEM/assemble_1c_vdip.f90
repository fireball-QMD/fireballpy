! ===========================================================================
! This routine assembles all of the one-center exchange-correlation
! interactions. The results are stored in vxc_1c and etotxc_1c.
! ===========================================================================
subroutine assemble_1c_vdip (iforce)
  use M_fdata
  implicit none
 
  integer, intent(in) :: iforce
 
  integer iatom
  integer imu, inu, ialpha, ibeta
  integer in1
  integer matom
  integer ixc
  integer indx

  real dccexc_1c
  real exc_1c
  real muexc_1c
  real I
  real Integral
  real, dimension (numorb_max,numorb_max) :: mu1xc
 
  Vdip_1c = 0.0d0
  dc_v_intra_dip_1c = 0.0d0
  do iatom = 1,natoms
    in1 = imass(iatom)
    matom = neigh_self(iatom)
    do indx = 1,Nlines_vdip1c(in1)
      !Here imu,inu,ialpha,ibeta are assumed to be pairwise different
      imu    = muR(indx,in1)
      inu    = nuR(indx,in1)
      ialpha = alphaR(indx,in1)
      ibeta  = betaR(indx,in1)
      I      = IR(indx,in1)
      if (imu .ne. inu) then
      !eq2?
        Vdip_1c(imu,inu,iatom) = (rho(ialpha,ibeta,matom,iatom)+rho(ibeta,ialpha,matom,iatom))*I*eq2
        Vdip_1c(inu,imu,iatom) = (rho(ialpha,ibeta,matom,iatom)+rho(ibeta,ialpha,matom,iatom))*I*eq2
        dc_v_intra_dip_1c = dc_v_intra_dip_1c - 0.5*(rho(imu,inu,matom,iatom)*Vdip_1c(imu,inu,iatom)+ rho(inu,imu,matom,iatom)*Vdip_1c(inu,imu,iatom))
      end if
      !There's another DC term remaining..!!  --> -xi^2 xi xi arising from rho_0*rho_dip   
    end do !end do indx
  end do !end do iatom

return
end subroutine assemble_1c_vdip
 
