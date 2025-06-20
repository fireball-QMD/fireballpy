subroutine Dassemble_ca_3c () 
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: smt_elect, natoms, ratom, imass, dip, dipp, neigh_comb, neigh_comj, neigh_comm, neigh_comn, numorb_max, sp_mat, &
    & Qin, rho, s_mat, xl, f3caa, f3cab, f3cac
  use M_fdata, only: nssh,rcutoff,Qneutral,lssh,num_orb
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer icount
  integer icount_sav
  integer imu
  integer in1
  integer in2
  integer indna
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer issh
  integer ix
  integer jatom
  integer jbeta
  integer jcount
  integer jcount_sav
  integer jssh
  integer mneigh
  real(double) cost
  real(double) distance13
  real(double) distance23
  real(double) dq1
  real(double) dq2
  real(double) dq3
  real(double) dstn_temp1
  real(double) dstn_temp2
  real(double) dterm
  real(double) dxn
  real(double) rcutoff_ialp
  real(double) rend1
  real(double) rend2
  real(double) sterm
  real(double) stn_temp1
  real(double) stn_temp2
  real(double) x
  real(double) y
  real(double) rcutoff_i
  real(double) rcutoff_j
  real(double), dimension (numorb_max, numorb_max) :: bcca
  real(double), dimension (numorb_max, numorb_max) :: bccax
  real(double), dimension (3, numorb_max, numorb_max) :: demnplA
  real(double), dimension (3, numorb_max, numorb_max) :: demnplB
  real(double), dimension (3, numorb_max, numorb_max) :: demnplC
  real(double), dimension (3, 3, 3) :: depsA
  real(double), dimension (3, 3, 3) :: depsB
  real(double), dimension (3) :: dpterm
  real(double), dimension (numorb_max, numorb_max) :: dstn1
  real(double), dimension (numorb_max, numorb_max) :: dstn2
  real(double), dimension (3, numorb_max, numorb_max) :: dstnA
  real(double), dimension (3, numorb_max, numorb_max) :: dstnB
  real(double), dimension (3, numorb_max, numorb_max) :: dstnC
  real(double), dimension (numorb_max, numorb_max) :: emnpl
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXa
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXb
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXc
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXa_sorp
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXb_sorp
  real(double), dimension (3, numorb_max, numorb_max) :: f3caXc_sorp
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r21
  real(double), dimension (3) :: rhat
  real(double), dimension (3) :: rhatA1
  real(double), dimension (3) :: rhatA2
  real(double), dimension (3) :: rna
  real(double), dimension (3) :: rnabc
  real(double), dimension (3) :: sighat
  real(double), dimension (3) :: spterm
  real(double), dimension (numorb_max, numorb_max) :: stn1
  real(double), dimension (numorb_max, numorb_max) :: stn2
  f3caa = 0.0d0
  f3cab = 0.0d0
  f3cac = 0.0d0
  do ialp = 1, natoms
   rna(:) = ratom(:,ialp)
   indna = imass(ialp)
   rcutoff_ialp = 0.0d0
   do imu = 1, nssh(indna)
    if (rcutoff(indna,imu) .gt. rcutoff_ialp)  rcutoff_ialp = rcutoff(indna,imu)
   end do
   do ineigh = 1, neigh_comn(ialp)
    mneigh = neigh_comm(ineigh,ialp)
    if (mneigh .ne. 0) then
     iatom = neigh_comj(1,ineigh,ialp)
     ibeta = neigh_comb(1,ineigh,ialp)
     r1(:) = ratom(:,iatom) + xl(:,ibeta)
     in1 = imass(iatom)
     dq1 = 0.0d0
     do issh = 1, nssh(in1)
      dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
     end do
    rcutoff_i = 0
    do imu = 1, nssh(in1)
     if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do
     jatom = neigh_comj(2,ineigh,ialp)
     jbeta = neigh_comb(2,ineigh,ialp)
     r2(:) = ratom(:,jatom) + xl(:,jbeta)
     in2 = imass(jatom)
     dq2 = 0.0d0
     do issh = 1, nssh(in2)
      dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
     end do
    rcutoff_j = 0
    do imu = 1, nssh(in2)
     if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
    end do
     r21(:) = r2(:) - r1(:)
     y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
     if (y .lt. 1.0d-05) then
      sighat(1) = 0.0d0
      sighat(2) = 0.0d0
      sighat(3) = 1.0d0
     else
      sighat(:) = r21(:)/y
     end if
     rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
     x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
     if (x .lt. 1.0d-03) then
      rhat(1) = 0.0d0
      rhat(2) = 0.0d0
      rhat(3) = 0.0d0
     else
      rhat(:) = rnabc(:)/x
     end if
     cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
     call epsilon (rhat, sighat, eps)
     call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
     distance13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2 + (rna(3) - r1(3))**2)
     distance23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2 + (rna(3) - r2(3))**2)
     if (distance13 .lt. 1.0d-05) distance13 = 1.0d-5
     if (distance23 .lt. 1.0d-05) distance23 = 1.0d-5
    rhatA1(:) = (rna(:) - r1(:))/distance13
    rhatA2(:) = (rna(:) - r2(:))/distance23
     dq3 = 0.0d0
     do issh = 1, nssh(indna)
      dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
     end do
     icount_sav = 0
     do issh = 1, nssh(in1)
      jcount_sav = 0
      do jssh = 1, nssh(in2)
       rend1 = rcutoff_i + rcutoff_ialp
       rend2 = rcutoff_j + rcutoff_ialp
       call smoother (distance13, rend1, smt_elect, stn_temp1, dstn_temp1)
       call smoother (distance23, rend2, smt_elect, stn_temp2, dstn_temp2)
       stn_temp1 = stn_temp1*stn_temp2
       stn_temp2 = 1.0d0 - stn_temp1
       do inu = 1, lssh(issh,in1)*2 + 1
        icount = icount_sav + inu
        do imu = 1, lssh(jssh,in2)*2 + 1
         jcount = jcount_sav + imu
         stn1(icount,jcount) = stn_temp1
         stn2(icount,jcount) = stn_temp2
         dstn1(icount,jcount) = dstn_temp1
         dstn2(icount,jcount) = dstn_temp2
        end do
       end do
       jcount_sav = jcount
      end do
      icount_sav = icount
     end do
     dstnB = 0.0d0
     dstnC = 0.0d0
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       dstnB(:,imu,inu) = - dstn1(imu,inu)*stn2(imu,inu)*rhatA1(:)
       dstnC(:,imu,inu) = - stn1(imu,inu)*dstn2(imu,inu)*rhatA2(:)
      end do
     end do
     dstnA = - dstnB - dstnC
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
            !if (y .lt. 1d-8) cycle  ! No sé que mierdas pasa
       sterm = dq3*s_mat(imu,inu,mneigh,iatom)/2.0d0
       dterm = dq3*dip(imu,inu,mneigh,iatom)/y
       spterm(:) = dq3*sp_mat(:,imu,inu,mneigh,iatom)/2.0d0
       dpterm(:) = dq3*(dipp(:,imu,inu,mneigh,iatom)/y + dip(imu,inu,mneigh,iatom)*sighat(:)/(y*y))
       emnpl(imu,inu) = (sterm - dterm)/distance13  + (sterm + dterm)/distance23
       demnplA(:,imu,inu) = - (sterm - dterm)*rhatA1(:)/distance13**2 - (sterm + dterm)*rhatA2(:)/distance23**2
       demnplB(:,imu,inu) = (sterm - dterm)*(rhatA1(:)/distance13**2) + (spterm(:) - dpterm(:))/distance13   + (spterm(:) + dpterm(:))/distance23 
       demnplC(:,imu,inu) = - demnplA(:,imu,inu) - demnplB(:,imu,inu)
      end do
     end do
     bcca = 0.0d0
     f3caXa = 0.0d0
     f3caXb = 0.0d0
     f3caXc = 0.0d0
     interaction = 1
     do isorp = 1, nssh(indna)
      call Dtrescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bccax, f3caXa_sorp, f3caXb_sorp, f3caXc_sorp)
      dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
      do inu = 1, num_orb(in2)
       do imu = 1, num_orb(in1)
        bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
        f3caXa(:,imu,inu) = f3caXa(:,imu,inu) + f3caXa_sorp(:,imu,inu)*dxn
        f3caXb(:,imu,inu) = f3caXb(:,imu,inu) + f3caXb_sorp(:,imu,inu)*dxn
        f3caXc(:,imu,inu) = f3caXc(:,imu,inu) + f3caXc_sorp(:,imu,inu)*dxn
       end do
      end do
     end do
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        f3caa(ix,ialp)  = f3caa(ix,ialp) + 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1(imu,inu)*f3caXa(ix,imu,inu) - dstnA(ix,imu,inu)*bcca(imu,inu) - stn2(imu,inu)*demnplA(ix,imu,inu)  + dstnA(ix,imu,inu)*emnpl(imu,inu))
        f3cab(ix,iatom) = f3cab(ix,iatom)+ 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1(imu,inu)*f3caXb(ix,imu,inu) - dstnB(ix,imu,inu)*bcca(imu,inu) - stn2(imu,inu)*demnplB(ix,imu,inu)  + dstnB(ix,imu,inu)*emnpl(imu,inu))
        f3cac(ix,jatom) = f3cac(ix,jatom)+ 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1(imu,inu)*f3caXc(ix,imu,inu) - dstnC(ix,imu,inu)*bcca(imu,inu) - stn2(imu,inu)*demnplC(ix,imu,inu)  + dstnC(ix,imu,inu)*emnpl(imu,inu))   
       end do
      end do
     end do   
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        f3caa(ix,ialp) = f3caa(ix,ialp)   + 2*rho(imu,inu,mneigh,iatom)*demnplA(ix,imu,inu)*eq2
        f3cab(ix,iatom) = f3cab(ix,iatom) + 2*rho(imu,inu,mneigh,iatom)*demnplB(ix,imu,inu)*eq2
        f3cac(ix,jatom) = f3cac(ix,jatom) + 2*rho(imu,inu,mneigh,iatom)*demnplC(ix,imu,inu)*eq2
       end do ! do ix
      end do ! do imu
     end do ! do inu
    end if
   end do
  end do
  return
end subroutine Dassemble_ca_3c
