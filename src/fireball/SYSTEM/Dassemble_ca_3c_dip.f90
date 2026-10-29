subroutine Dassemble_ca_3c_dip () 
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: eq2
  use M_system, only: smt_elect, natoms, ratom, imass, neigh_comb, neigh_comj, neigh_comm, neigh_comn, numorb_max, sp_mat, &
    & dippc, Qin, rho, s_mat, dipc, xl, f3caa, f3cab, f3cac
  use M_fdata, only: nssh,rcutoff,Qneutral,num_orb,nspecies
  implicit none
  integer ialp
  integer iatom
  integer ibeta
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
  integer mneigh
  real(dp) cost
  real(dp) distance13
  real(dp) distance23
  real(dp) dq1
  real(dp) dq2
  real(dp) dq3
  real(dp) dterm
  real(dp) dxn
  real(dp) rcutoff_ialp
  real(dp) rend1
  real(dp) rend2
  real(dp) sterm
  real(dp) stn_temp1
  real(dp) stn_temp2
  real(dp) x
  real(dp) y
  real(dp) rcutoff_i
  real(dp) rcutoff_j
  real(dp), dimension (numorb_max, numorb_max) :: bcca
  real(dp), dimension (numorb_max, numorb_max) :: bccax
  real(dp), dimension (3, numorb_max, numorb_max) :: demnplA
  real(dp), dimension (3, numorb_max, numorb_max) :: demnplB
  real(dp), dimension (3, numorb_max, numorb_max) :: demnplC
  real(dp), dimension (3, 3, 3) :: depsA
  real(dp), dimension (3, 3, 3) :: depsB
  real(dp) dstn1
  real(dp) dstn2
  real(dp), dimension (3) :: dstnA
  real(dp), dimension (3) :: dstnB
  real(dp), dimension (3) :: dstnC
  real(dp), dimension (numorb_max, numorb_max) :: emnpl
  real(dp), dimension (3, 3) :: eps
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXa
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXb
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXc
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXa_sorp
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXb_sorp
  real(dp), dimension (3, numorb_max, numorb_max) :: f3caXc_sorp
  real(dp), dimension (3) :: r1
  real(dp), dimension (3) :: r2
  real(dp), dimension (3) :: r21
  real(dp), dimension (3) :: rhat
  real(dp), dimension (3) :: rhatA1
  real(dp), dimension (3) :: rhatA2
  real(dp), dimension (3) :: rna
  real(dp), dimension (3) :: rnabc
  real(dp), dimension (3) :: sighat
  real(dp), dimension (3) :: spterm
  real(dp), dimension (3) :: ddterm
  real(dp), dimension (3) :: dptermA
  real(dp), dimension (3) :: dptermB
  real(dp) stn1
  real(dp) stn2
  f3caa = 0.0_dp
  f3cab = 0.0_dp
  f3cac = 0.0_dp
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    rcutoff_ialp = 0.0_dp
    do imu = 1, nssh(indna)
      if (rcutoff(imu,indna) .gt. rcutoff_ialp) rcutoff_ialp = rcutoff(imu,indna)
    end do
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh,ialp)
      if (mneigh .ne. 0) then
        iatom = neigh_comj(1,ineigh,ialp)
        ibeta = neigh_comb(1,ineigh,ialp)
        r1(:) = ratom(:,iatom) + xl(:,ibeta)
        in1 = imass(iatom)
        dq1 = 0.0_dp
        do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
        end do
        rcutoff_i = 0
        do imu = 1, nssh(in1)
          if (rcutoff(imu,in1) .gt. rcutoff_i) rcutoff_i = rcutoff(imu,in1)
        end do
        jatom = neigh_comj(2,ineigh,ialp)
        jbeta = neigh_comb(2,ineigh,ialp)
        r2(:) = ratom(:,jatom) + xl(:,jbeta)
        in2 = imass(jatom)
        dq2 = 0.0_dp
        do issh = 1, nssh(in2)
          dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
        end do
        rcutoff_j = 0
        do imu = 1, nssh(in2)
          if (rcutoff(imu,in2) .gt. rcutoff_j) rcutoff_j = rcutoff(imu,in2)
        end do
        r21(:) = r2(:) - r1(:)
        y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
        if (y .lt. 1.0e-05_dp) then
          sighat(1) = 0.0_dp
          sighat(2) = 0.0_dp
          sighat(3) = 1.0_dp
        else
          sighat(:) = r21(:)/y
        end if
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0_dp)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        if (x .lt. 1.0e-03_dp) then
          rhat(1) = 0.0_dp
          rhat(2) = 0.0_dp
          rhat(3) = 0.0_dp
        else
          rhat(:) = rnabc(:)/x
        end if
        cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
        call epsilon (rhat, sighat, eps)
        call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
        distance13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2  + (rna(3) - r1(3))**2)
        distance23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2  + (rna(3) - r2(3))**2)
        if (distance13 .lt. 1.0e-05_dp) distance13 = 1.0e-5_dp
        if (distance23 .lt. 1.0e-05_dp) distance23 = 1.0e-5_dp
        rhatA1(:) = (rna(:) - r1(:))/distance13
        rhatA2(:) = (rna(:) - r2(:))/distance23
        dq3 = 0.0_dp
        do issh = 1, nssh(indna)
           dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
        end do
        rend1 = rcutoff_i + rcutoff_ialp
        rend2 = rcutoff_j + rcutoff_ialp
        call smoother (distance13, rend1, smt_elect, stn_temp1, dstn1)
        call smoother (distance23, rend2, smt_elect, stn_temp2, dstn2)
        stn1 = stn_temp1*stn_temp2
        stn2 = 1.0_dp - stn1
        dstnB = 0.0_dp
        dstnC = 0.0_dp
        dstnB(:) = - dstn1*stn2*rhatA1(:)
        dstnC(:) = - stn1*dstn2*rhatA2(:)
        dstnA = - dstnB - dstnC
        if (x .lt. 1.0e-05_dp) then
          emnpl = 0.0_dp
          demnplA = 0.0_dp
          demnplB = 0.0_dp
          demnplC = 0.0_dp
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,mneigh,iatom)
              dterm = (dipc(1,imu,inu,mneigh,iatom)*rnabc(1)  + dipc(2,imu,inu,mneigh,iatom)*rnabc(2)   + dipc(3,imu,inu,mneigh,iatom)*rnabc(3))
              ddterm(:) = (dippc(:,1,imu,inu,mneigh,iatom)*rnabc(1)  + dippc(:,2,imu,inu,mneigh,iatom)*rnabc(2)   + dippc(:,3,imu,inu,mneigh,iatom)*rnabc(3))
              spterm(:) = sp_mat(:,imu,inu,mneigh,iatom)
              dptermA(:)=   dipc(:,imu,inu,mneigh,iatom)/(x*x*x)  - 3*dterm*rnabc(:)/(x*x*x*x*x)
              dptermB(:)= - 0.50*dipc(:,imu,inu,mneigh,iatom)/(x*x*x)  + 0.50*3*dterm*rnabc(:)/(x*x*x*x*x)  + ddterm(:)/(x*x*x)
              emnpl(imu,inu) =  dq3*sterm/x  + dq3*dterm/(x*x*x)
              demnplA(:,imu,inu) = - dq3*sterm*rnabc(:)/(x*x*x)  + dq3*dptermA(:)
              demnplB(:,imu,inu) =  dq3*0.50*sterm*rnabc(:)/(x*x*x)  + dq3*spterm(:)/(x)  + dq3*dptermB(:)
              demnplC(:,imu,inu) = - demnplA(:,imu,inu) - demnplB(:,imu,inu)
            end do
          end do
        end if
        bcca = 0.0_dp
        f3caXa = 0.0_dp
        f3caXb = 0.0_dp
        f3caXc = 0.0_dp
        interaction = 1
        do isorp = 1, nssh(indna)
          call Dtrescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bccax, f3caXa_sorp, f3caXb_sorp, f3caXc_sorp, nspecies)
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
              f3caa(ix,ialp) = f3caa(ix,ialp)   + 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1*f3caXa(ix,imu,inu) - dstnA(ix)*bcca(imu,inu) - stn2*demnplA(ix,imu,inu)  + dstnA(ix)*emnpl(imu,inu))
              f3cab(ix,iatom) = f3cab(ix,iatom) + 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1*f3caXb(ix,imu,inu) - dstnB(ix)*bcca(imu,inu) - stn2*demnplB(ix,imu,inu)  + dstnB(ix)*emnpl(imu,inu))
              f3cac(ix,jatom) = f3cac(ix,jatom) + 2*rho(imu,inu,mneigh,iatom)*eq2*(stn1*f3caXc(ix,imu,inu) - dstnC(ix)*bcca(imu,inu) - stn2*demnplC(ix,imu,inu)  + dstnC(ix)*emnpl(imu,inu))   
            end do
          end do
        end do   
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            do ix = 1, 3
              f3caa(ix,ialp) = f3caa(ix,ialp)  + 2*rho(imu,inu,mneigh,iatom)*demnplA(ix,imu,inu)*eq2
              f3cab(ix,iatom) = f3cab(ix,iatom) + 2*rho(imu,inu,mneigh,iatom)*demnplB(ix,imu,inu)*eq2
              f3cac(ix,jatom) = f3cac(ix,jatom) + 2*rho(imu,inu,mneigh,iatom)*demnplC(ix,imu,inu)*eq2
            end do
          end do 
        end do 
      end if !mneigh .ne. 0
    end do !ineigh
  end do !atoms   
end subroutine Dassemble_ca_3c_dip
