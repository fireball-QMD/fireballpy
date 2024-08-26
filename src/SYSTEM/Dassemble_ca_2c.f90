subroutine Dassemble_ca_2c ()
  use M_constants, only: wp, eq2
 use M_system
 use M_fdata, only: nssh, Qneutral, rcutoff, lssh, num_orb
 implicit none
 integer iatom
 integer icount
 integer icount_sav
 integer ierror
 integer imu
 integer in1
 integer in2
 integer in3
 integer ineigh
 integer interaction
 integer inu
 integer isorp
 integer issh
 integer ix
 integer jatom
 integer jcount
 integer jcount_sav
 integer jssh
 integer kforce
 integer matom
 integer mbeta
 real(wp) dq1
 real(wp) dq2
 real(wp) dstn_temp
 real(wp) dterm_1
 real(wp) dterm_2
 real(wp) dxn
 real(wp) rcutoff_j
 real(wp) rend
 real(wp) rend1
 real(wp) rend2
 real(wp) sterm_1
 real(wp) sterm_2
 real(wp) stn_temp1
 real(wp) stn_temp2
 real(wp) y
 real(wp) rcutoff_i
 real(wp), dimension (numorb_max, numorb_max) :: bcca
 real(wp), dimension (3, numorb_max, numorb_max) :: bccap
 real(wp), dimension (3, numorb_max, numorb_max) :: bccapx
 real(wp), dimension (numorb_max, numorb_max) :: bccax
 real(wp), dimension (numorb_max, numorb_max) :: demnpl
 real(wp), dimension (3, 3, 3) :: deps
 real(wp), dimension (3, numorb_max, numorb_max) :: dewaldsr
 real(wp), dimension (3) :: dpterm_1
 real(wp), dimension (3) :: dpterm_2
 real(wp), dimension (numorb_max, numorb_max) :: dstn1
 real(wp), dimension (numorb_max, numorb_max) :: dstn2
 real(wp), dimension (numorb_max, numorb_max) :: emnpl
 real(wp), dimension (3, 3) :: eps
 real(wp), dimension (3) :: r1
 real(wp), dimension (3) :: r2
 real(wp), dimension (3) :: r21
 real(wp), dimension (3) :: sighat
 real(wp), dimension (3) :: spterm_1
 real(wp), dimension (3) :: spterm_2
 real(wp), dimension (numorb_max, numorb_max) :: stn1
 real(wp), dimension (numorb_max, numorb_max) :: stn2
 faca = 0.0d0
 fotca = 0.0d0
 do iatom = 1, natoms
  matom = neigh_self(iatom)
  r1(:) = ratom(:,iatom)
  in1 = imass(iatom)
   rcutoff_i = 0
   do imu = 1, nssh(in1)
    if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
   end do
  dq1 = 0.0
  do issh = 1, nssh(in1)
   dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
  end do
  do ineigh = 1, neighn(iatom)
   mbeta = neigh_b(ineigh,iatom)
   jatom = neigh_j(ineigh,iatom)
   r2(:) = ratom(:,jatom) + xl(:,mbeta)
   in2 = imass(jatom)
   rcutoff_j = 0
   do imu = 1, nssh(in2)
    if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
   end do
   dq2 = 0.0d0
   do issh = 1, nssh(in2)
    dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
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
   call epsilon (r2, sighat, eps)
   call deps2cent (r1, r2, eps, deps)
   stn1 = 1.0d0
   dstn1 = 0.0d0
   stn2 = 0.0d0
   dstn2 = 0.0d0
   emnpl = 0.0d0
   demnpl = 0.0d0
   dewaldsr = 0.0d0
   if (y .gt. 1.0d-4) then
    icount_sav = 0
    do issh = 1, nssh(in1)
     jcount_sav = 0
     do jssh = 1, nssh(in1)
      rend = rcutoff_i + rcutoff_j
      call smoother (y, rend, smt_elect, stn_temp1, dstn_temp)
      stn_temp2 = 1.0d0 - stn_temp1
      do inu = 1, lssh(issh,in1)*2 + 1
       icount = icount_sav + inu
       do imu = 1, lssh(jssh,in1)*2 + 1
        jcount = jcount_sav + imu
        stn1(icount,jcount) = stn_temp1
        stn2(icount,jcount) = stn_temp2
        dstn1(icount,jcount) = dstn_temp
        dstn2(icount,jcount) = - dstn_temp
       end do
      end do
      jcount_sav = jcount
     end do
     icount_sav = icount
    end do
    do inu = 1, num_orb(in1)
     do imu = 1, num_orb(in1)
      emnpl(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)*dq2
      demnpl(imu,inu) = - (s_mat(imu,inu,matom,iatom)/(y*y))*dq2
      dewaldsr(:,imu,inu) = - demnpl(imu,inu)*sighat(:)*eq2
     end do
    end do
   end if
   bcca = 0.0d0
   bccap = 0.0d0
   kforce = 1
   interaction = 4
   in3 = in1
   do isorp = 1, nssh(in2)
    call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  eps, deps, bccax, bccapx)
    dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
    do inu = 1, num_orb(in3)
     do imu = 1, num_orb(in1)
      bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
      bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
     end do
    end do
   end do
   do inu = 1, num_orb(in3)
    do imu = 1, num_orb(in1)
     bccap(:,imu,inu) =  stn1(imu,inu)*bccap(:,imu,inu) - dstn1(imu,inu)*bcca(imu,inu)*sighat(:)  - (stn2(imu,inu)*demnpl(imu,inu)  + dstn2(imu,inu)*emnpl(imu,inu))*sighat(:)
    end do
   end do
   do inu = 1, num_orb(in3)
    do imu = 1, num_orb(in1)
     do ix = 1, 3
      faca(ix,ineigh,iatom) = faca(ix,ineigh,iatom)  - rho(imu,inu,matom,iatom)*bccap(ix,imu,inu)*eq2  + rho(imu,inu,matom,iatom)*dewaldsr(ix,imu,inu)
     end do
    end do
   end do
   if (iatom .eq. jatom .and. mbeta .eq. 0) then
   else
    do inu = 1, num_orb(in2)
     do imu = 1, num_orb(in1)
      sterm_1 = dq1*eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
      sterm_2 = dq2*eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
      dterm_1 = dq1*eq2*dip(imu,inu,ineigh,iatom)/(y*y)
      dterm_2 = dq2*eq2*dip(imu,inu,ineigh,iatom)/(y*y)
      spterm_1(:) = dq1*eq2*sp_mat(:,imu,inu,ineigh,iatom)/(2.0d0*y)
      spterm_2(:) = dq2*eq2*sp_mat(:,imu,inu,ineigh,iatom)/(2.0d0*y)
      dpterm_1(:) = dq1*eq2*dipp(:,imu,inu,ineigh,iatom)/(y*y)
      dpterm_2(:) = dq2*eq2*dipp(:,imu,inu,ineigh,iatom)/(y*y)
      dewaldsr(:,imu,inu) =  + (spterm_1(:) + dpterm_1(:)) + (sterm_1 + 2.0d0*dterm_1)*sighat(:)/y  + (spterm_2(:) - dpterm_2(:)) + (sterm_2 - 2.0d0*dterm_2)*sighat(:)/y
     end do
    end do
    bccap = 0.0d0
    interaction = 2
    in3 = in2
    do isorp = 1, nssh(in1)
     call doscentros (interaction, isorp, kforce, in1, in1, in3, y,eps, deps, bccax, bccapx)
     dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
     do imu = 1, num_orb(in1)
      do inu = 1, num_orb(in3)
       bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
      end do
     end do
    end do
    do inu = 1, num_orb(in3)
     do imu = 1, num_orb(in1)
      do ix = 1, 3
       fotca(ix,ineigh,iatom) = fotca(ix,ineigh,iatom)  - rho(imu,inu,ineigh,iatom)*bccap(ix,imu,inu)*eq2  + 0.5d0*rho(imu,inu,ineigh,iatom)*dewaldsr(ix,imu,inu)
      end do
     end do
    end do
   end if
  end do
 end do
 return
end

