! This routine assembles all of the two-center and degenerate two-center interactions for DOGS.
subroutine assemble_ca_2c ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system, only: iforce, smt_elect, natoms, ratom, imass, ewaldsr, dip, dipp, neigh_b, neigh_j, neighn, neigh_self, &
    & numorb_max, Qin, s_mat, vca, xl
  use M_fdata, only: nssh, num_orb, rcutoff, Qneutral, lssh
  implicit none
  integer(c_long) iatom
  integer(c_long) icount
  integer(c_long) icount_sav
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) in3
  integer(c_long) ineigh
  integer(c_long) interaction
  integer(c_long) inu
  integer(c_long) isorp
  integer(c_long) issh
  integer(c_long) jatom
  integer(c_long) jcount
  integer(c_long) jcount_sav
  integer(c_long) jssh
  integer(c_long) kforce
  integer(c_long) matom
  integer(c_long) mbeta
  real(c_double) dq1
  real(c_double) dq2
  real(c_double) dterm_1
  real(c_double) dterm_2
  real(c_double) dstn_temp
  real(c_double) dxn
  real(c_double) rcutoff_j
  real(c_double) rend
  real(c_double) sterm_1
  real(c_double) sterm_2
  real(c_double) stn_temp1
  real(c_double) stn_temp2
  real(c_double) y
  real(c_double) rcutoff_i
  real(c_double), dimension (numorb_max, numorb_max) :: bcca
  real(c_double), dimension (3, numorb_max, numorb_max) :: bccapx
  real(c_double), dimension (numorb_max, numorb_max) :: bccax
  real(c_double), dimension (3, 3, 3) :: deps
  real(c_double), dimension (numorb_max, numorb_max) :: dipx
  real(c_double), dimension (3, numorb_max, numorb_max) :: dippx
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl_noq
  real(c_double), dimension (3, 3) :: eps
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r21
  real(c_double), dimension (3) :: sighat
  real(c_double), dimension (numorb_max, numorb_max) :: stn1
  real(c_double), dimension (numorb_max, numorb_max) :: stn2

  vca = 0.0d0
  ewaldsr = 0.0d0
  dip = 0.0d0
  dipp = 0.0d0

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
 
      ! CALL DOSCENTROS AND GET DIP
      isorp = 0
      interaction = 9
      in3 = in2
      call doscentros (interaction, isorp, iforce, in1, in2, in3, y,  eps, deps, dipx, dippx)
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          dip(imu,inu,ineigh,iatom) = dipx(imu,inu)
          if (iforce .eq. 1) dipp(:,imu,inu,ineigh,iatom) = dippx(:,imu,inu)
        end do
      end do
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
