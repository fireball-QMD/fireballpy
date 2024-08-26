subroutine Dassemble_ca_3c_dip () 
  use M_constants, only: wp, eq2
  use M_system
  use M_fdata, only: isorpmax, nssh,rcutoff,Qneutral,num_orb,nspecies
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer icount
  integer ierror
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
  real(wp) cost
  real(wp) distance13
  real(wp) distance23
  real(wp) dq1
  real(wp) dq2
  real(wp) dq3
  real(wp) dstn_temp1
  real(wp) dstn_temp2
  real(wp) dterm
  real(wp) dxn
  real(wp) rcutoff_ialp
  real(wp) rend1
  real(wp) rend2
  real(wp) sterm
  real(wp) stn_temp1
  real(wp) stn_temp2
  real(wp) x
  real(wp) y
  real(wp) rcutoff_i
  real(wp) rcutoff_j
  real(wp), dimension (numorb_max, numorb_max) :: bcca
  real(wp), dimension (numorb_max, numorb_max) :: bccax
  real(wp), dimension (3, numorb_max, numorb_max) :: demnplA
  real(wp), dimension (3, numorb_max, numorb_max) :: demnplB
  real(wp), dimension (3, numorb_max, numorb_max) :: demnplC
  real(wp), dimension (3, 3, 3) :: depsA
  real(wp), dimension (3, 3, 3) :: depsB
  real(wp), dimension (3) :: dpterm
  real(wp) dstn1
  real(wp) dstn2
  real(wp), dimension (3) :: dstnA
  real(wp), dimension (3) :: dstnB
  real(wp), dimension (3) :: dstnC
  real(wp), dimension (numorb_max, numorb_max) :: emnpl
  real(wp), dimension (3, 3) :: eps
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXa
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXb
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXc
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXa_sorp
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXb_sorp
  real(wp), dimension (3, numorb_max, numorb_max) :: f3caXc_sorp
  real(wp), dimension (3) :: r1
  real(wp), dimension (3) :: r2
  real(wp), dimension (3) :: r21
  real(wp), dimension (3) :: rhat
  real(wp), dimension (3) :: rhatA1
  real(wp), dimension (3) :: rhatA2
  real(wp), dimension (3) :: rna
  real(wp), dimension (3) :: rnabc
  real(wp), dimension (3) :: sighat
  real(wp), dimension (3) :: spterm
  real(wp), dimension (3) :: ddterm
  real(wp), dimension (3) :: dptermA
  real(wp), dimension (3) :: dptermB
  real(wp) stn1
  real(wp) stn2
  f3caa = 0.0d0
  f3cab = 0.0d0
  f3cac = 0.0d0
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    rcutoff_ialp = 0.0d0
    do imu = 1, nssh(indna)
      if (rcutoff(indna,imu) .gt. rcutoff_ialp) rcutoff_ialp = rcutoff(indna,imu)
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
          write (*,*) ' There is an error here in assemble_3c.f '
          write (*,*) ' r1 = r2!!!! BAD!!!! '
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
        distance13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2  + (rna(3) - r1(3))**2)
        distance23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2  + (rna(3) - r2(3))**2)
        if (distance13 .gt. 1.0d-05) then
          rhatA1(:) = (rna(:) - r1(:))/distance13
        else
          write (*,*) ' distance13 is too small in Dassemble_ca_3c_dip.f '
          write (*,*) ' This can not be so!!!! '
        end if
        if (distance23 .gt. 1.0d-05) then
            rhatA2(:) = (rna(:) - r2(:))/distance23
        else
          write (*,*) ' distance23 is too small in Dassemble_ca_3c_dip.f '
          write (*,*) ' This can not be so!!!! '
        end if
        dq3 = 0.0d0
        do issh = 1, nssh(indna)
           dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
        end do
        rend1 = rcutoff_i + rcutoff_ialp
        rend2 = rcutoff_j + rcutoff_ialp
        call smoother (distance13, rend1, smt_elect, stn_temp1, dstn1)
        call smoother (distance23, rend2, smt_elect, stn_temp2, dstn2)
        stn1 = stn_temp1*stn_temp2
        stn2 = 1.0d0 - stn1
        dstnB = 0.0d0
        dstnC = 0.0d0
        dstnB(:) = - dstn1*stn2*rhatA1(:)
        dstnC(:) = - stn1*dstn2*rhatA2(:)
        dstnA = - dstnB - dstnC
        if (x .lt. 1.0d-05) then
          emnpl = 0.0d0
          demnplA = 0.0d0
          demnplB = 0.0d0
          demnplC = 0.0d0
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
        bcca = 0.0d0
        f3caXa = 0.0d0
        f3caXb = 0.0d0
        f3caXc = 0.0d0
        interaction = 1
        do isorp = 1, nssh(indna)
          call Dtrescentros (interaction, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bccax, f3caXa_sorp, f3caXb_sorp, f3caXc_sorp, nspecies)
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
end
  
