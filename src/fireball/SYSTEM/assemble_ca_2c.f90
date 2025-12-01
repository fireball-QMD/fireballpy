! This routine assembles all of the two-center and degenerate two-center interactions for DOGS.
subroutine assemble_ca_2c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: smt_elect, natoms, ratom, imass, ewaldsr, dip, neigh_b, neigh_j, neighn, neigh_self, &
    & numorb_max, Qin, s_mat, vca, xl, g_h, kscf, iqout
  use M_fdata, only: nssh, num_orb, rcutoff, Qneutral, lssh
  implicit none
  integer iatom
  integer icount
  integer icount_sav
  integer imu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer issh
  integer jatom
  integer jcount
  integer jcount_sav
  integer jssh
  integer kforce
  integer matom
  integer mbeta
  integer alpha, beta
  real(double) aux
  real(double) dq1
  real(double) dq2
  real(double) dterm_1
  real(double) dterm_2
  real(double) dstn_temp
  real(double) dxn
  real(double) rcutoff_j
  real(double) rend
  real(double) sterm_1
  real(double) sterm_2
  real(double) stn_temp1
  real(double) stn_temp2
  real(double) y
  real(double) rcutoff_i
  real(double), dimension (numorb_max, numorb_max) :: bcca
  real(double), dimension (3, numorb_max, numorb_max) :: bccapx
  real(double), dimension (numorb_max, numorb_max) :: bccax
  real(double), dimension (3, 3, 3) :: deps
  real(double), dimension (numorb_max, numorb_max) :: emnpl
  real(double), dimension (numorb_max, numorb_max) :: emnpl_noq
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r21
  real(double), dimension (3) :: sighat
  real(double), dimension (numorb_max, numorb_max) :: stn1
  real(double), dimension (numorb_max, numorb_max) :: stn2

  vca = 0.0d0
  ewaldsr = 0.0d0
  if (Kscf .eq. 1 .and. iqout .eq. 6) then
    g_h  = 0.0d0
  end if
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    rcutoff_i = 0
    do imu = 1, nssh(in1)
      if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do
    dq1 = 0.0d0
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
 
      dq2 = 0.0d0
      do issh = 1, nssh(in2)
        dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      stn1 = 1.0d0
      stn2 = 0.0d0
      emnpl = 0.0d0
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
              end do
            end do
            jcount_sav = jcount
          end do
          icount_sav = icount
        end do
        do inu = 1, num_orb(in1)
          do imu = 1, num_orb(in1)
            emnpl(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)*dq2
            emnpl_noq(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)
            ewaldsr(imu,inu,matom,iatom) =  ewaldsr(imu,inu,matom,iatom) + emnpl(imu,inu)*eq2
            if (Kscf .eq. 1 .and. iqout .eq. 6) then
              do issh = 1, nssh(in2)
                g_h(imu,inu,issh,jatom,matom,iatom)  =  g_h(imu,inu,issh,jatom,matom,iatom) - emnpl_noq(imu,inu)*eq2
              end do 
            end if
          end do
        end do
      end if
      ! CALL DOSCENTROS AND GET VNA FOR ATOM CASE
      bcca = 0.0d0
      kforce = 0
      interaction = 4
      in3 = in1
      do isorp = 1, nssh(in2)
        call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bccax, bccapx)
        dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
            if (Kscf .eq. 1 .and. iqout .eq. 6) then
              g_h(imu,inu,isorp,jatom,matom,iatom)  =  g_h(imu,inu,isorp,jatom,matom,iatom) + (stn1(imu,inu)*bccax(imu,inu) + stn2(imu,inu)*emnpl_noq(imu,inu))*eq2
            end if
          end do
        end do
      end do
      do inu = 1, num_orb(in3)
        do imu = 1, num_orb(in1)
          vca(imu,inu,matom,iatom) = vca(imu,inu,matom,iatom) + (stn1(imu,inu)*bcca(imu,inu) + stn2(imu,inu)*emnpl(imu,inu))*eq2
        end do
      end do
      ! ASSEMBLE EWALDSR FOR ONTOP CASE
      if (y .gt. 1.0d-4) then
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            sterm_1 = eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
            sterm_2 = eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
            dterm_1 = eq2*dip(imu,inu,ineigh,iatom)/(y*y)
            dterm_2 = eq2*dip(imu,inu,ineigh,iatom)/(y*y)
            ewaldsr(imu,inu,ineigh,iatom) = ewaldsr(imu,inu,ineigh,iatom) + (dq1*sterm_1 + dq1*dterm_1) + (dq2*sterm_2 - dq2*dterm_2)
            if (Kscf .eq. 1 .and. iqout .eq. 6) then
              do issh = 1, nssh(in2)
                emnpl_noq(imu,inu) = sterm_1-dterm_1  ! on top right
                g_h(imu,inu,issh,jatom,ineigh,iatom)  =  g_h(imu,inu,issh,jatom,ineigh,iatom) - emnpl_noq(imu,inu)
              end do
              do issh = 1, nssh(in1)
                emnpl_noq(imu,inu) = sterm_1+dterm_1  ! on top left
                g_h(imu,inu,issh,jatom,ineigh,iatom)  =  g_h(imu,inu,issh,jatom,ineigh,iatom) - emnpl_noq(imu,inu)
              end do
            end if
          end do
        end do
      end if
   
      ! CALL DOSCENTROS AND GET VNA FOR ONTOP CASE
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        ! Do nothing here - special case. Interaction already calculated in atm case.
      else
        bcca = 0.0d0
        interaction = 2
        in3 = in2
        do isorp = 1, nssh(in1)
          call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, bccax, bccapx)
          dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
              if (Kscf .eq. 1 .and. iqout .eq. 6) then
                g_h(imu,inu,isorp,iatom,ineigh,iatom)  =  g_h(imu,inu,isorp,iatom,ineigh,iatom) + bccax(imu,inu)*eq2
              end if
            end do
          end do
        end do
        interaction = 3
        in3 = in2
        do isorp = 1, nssh(in2)
          call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bccax, bccapx)
          dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
              if (Kscf .eq. 1 .and. iqout .eq. 6) then
                g_h(imu,inu,isorp,jatom,ineigh,iatom)  =  g_h(imu,inu,isorp,jatom,ineigh,iatom) + bccax(imu,inu)*eq2
              end if
            end do
          end do
        end do
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            vca(imu,inu,ineigh,iatom) =  vca(imu,inu,ineigh,iatom) + bcca(imu,inu)*eq2
          end do
        end do
      end if
    end do   ! do ineigh
  end do   ! do iatom
  return
end subroutine assemble_ca_2c
