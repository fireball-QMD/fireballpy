! This routine assembles all of the two-center and degenerate two-center interactions for DOGS.
subroutine assemble_ca_2c_dip ()
  use M_system
  use M_fdata, only: nssh,rcutoff,Qneutral,num_orb
  use m_constants, only: eq2
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
  integer jatom
  integer jcount
  integer jcount_sav
  integer jssh
  integer kforce
  integer matom
  integer mbeta
  integer my_proc
  integer ix
  integer iy
 
  real dq1
  real dq2
  real dterm
  real dterm_1
  real dterm_2
  real dstn_temp
  real dxn
  real rcutoff_j
  real rend
  real rend1
  real rend2
  real sterm_1
  real sterm_2
  real y
  real rcutoff_i
 
  real, dimension (numorb_max, numorb_max) :: bcca
  real, dimension (3, numorb_max, numorb_max) :: bccapx
  real, dimension (numorb_max, numorb_max) :: bccax
  real, dimension (3, 3, 3) :: deps
  real, dimension (numorb_max, numorb_max) :: dipx
  real, dimension (3, numorb_max, numorb_max) :: dippx
  real, dimension (numorb_max, numorb_max) :: emnpl
  real, dimension (numorb_max, numorb_max) :: emnpl_noq
  real, dimension (3, 3) :: eps
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
  real, dimension (3) :: sighat
  real stn1
  real stn2
  vca = 0.0d0
  ewaldsr = 0.0d0
  if (Kscf .eq. 1 .and. iqout .eq. 6) then
    gvhxc = 0.0d0
  end if ! end if Kscf .eq. 1

  do iatom = 1,natoms
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
    do ineigh = 1, neighn(iatom)   ! <==== loop 2 over i's neighbors
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
      ! ASSEMBLE EWALDSR AND EMNPL FOR ATM CASE
      dq2 = 0.0d0
      do issh = 1, nssh(in2)
        dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      stn1 = 1.0d0
      stn2 = 0.0d0
      emnpl = 0.0d0
      emnpl_noq = 0.0d0
      if (y .gt. 1.0d-4) then
        rend = rcutoff_i + rcutoff_j
        call smoother (y, rend, smt_elect, stn1, dstn_temp)
        stn2 = 1.0d0 - stn1
        do inu = 1, num_orb(in1)
          do imu = 1, num_orb(in1)
            dterm = (dipc(1,imu,inu,matom,iatom)*r21(1)  + dipc(2,imu,inu,matom,iatom)*r21(2) + dipc(3,imu,inu,matom,iatom)*r21(3))
            emnpl(imu,inu) =  dq2*(s_mat(imu,inu,matom,iatom)/y) + dq2*(dterm/(y*y*y))
            emnpl_noq(imu,inu) =  (s_mat(imu,inu,matom,iatom)/y) + (dterm/(y*y*y))
            ewaldsr(imu,inu,matom,iatom) =  ewaldsr(imu,inu,matom,iatom) + emnpl(imu,inu)*eq2
            if (Kscf .eq. 1 .and. iqout .eq. 6) then
              do issh = 1, nssh(in2)
                gvhxc(imu,inu,issh,jatom,matom,iatom) =  gvhxc(imu,inu,issh,jatom,matom,iatom) - emnpl_noq(imu,inu)*eq2
              end do ! end do issh
            end if ! end if Kscf .eq. 1 .and. iqout .eq. 6
          end do
        end do
      end if
 
      ! CALL DOSCENTROS AND GET VNA FOR ATOM CASE
      bcca = 0.0d0
      kforce = 0
      interaction = 4
      in3 = in1
      do isorp = 1, nssh(in2)
        call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  eps, deps, bccax, bccapx)
        dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
            if (Kscf .eq. 1) then
              if (iqout .eq. 6) then
                gvhxc(imu,inu,isorp,jatom,matom,iatom) =  gvhxc(imu,inu,isorp,jatom,matom,iatom) +  (stn1*bccax(imu,inu) + stn2*emnpl_noq(imu,inu))*eq2
              end if ! end if iqout .eq. 6
            end if ! end if Kscf .eq. 1
          end do
        end do
      end do  ! isorp
      do inu = 1, num_orb(in3)
        do imu = 1, num_orb(in1)
          vca(imu,inu,matom,iatom) = vca(imu,inu,matom,iatom) + (stn1*bcca(imu,inu) + stn2*emnpl(imu,inu))*eq2
        end do
      end do
 
      !NOTE (true_dipoles, 2017):  We do not need to compute ewaldsr in the on-top case because we have excluded
      ! this case already in the new assemble_lr.f90 subroutine!!
      ! CALL DOSCENTROS AND GET VNA FOR ONTOP CASE
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        ! Do nothing here - special case. Interaction already calculated in atm case.
      else
        ! Initialize bcca for charged atom interactions.
        bcca = 0.0d0
        interaction = 2
        in3 = in2
        do isorp = 1, nssh(in1)
          call doscentros (interaction, isorp, kforce, in1, in1, in3, y,eps, deps, bccax, bccapx)
          dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
              if (Kscf .eq. 1) then
                if (iqout .eq. 6) then
                  gvhxc(imu,inu,isorp,iatom,ineigh,iatom) =gvhxc(imu,inu,isorp,iatom,ineigh,iatom) + bccax(imu,inu)*eq2
                end if ! end if iqout .eq. 6
              end if ! end if Kscf .eq. 1
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
              if (Kscf .eq. 1) then
                if (iqout .eq. 6) then
                  gvhxc(imu,inu,isorp,jatom,ineigh,iatom) =  gvhxc(imu,inu,isorp,jatom,ineigh,iatom) +  bccax(imu,inu)*eq2
                end if ! end if iqout .eq. 6
              end if ! end if Kscf .eq. 1
            end do
          end do
        end do
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
          vca(imu,inu,ineigh,iatom) = vca(imu,inu,ineigh,iatom) + bcca(imu,inu)*eq2
          end do
        end do
      end if
    end do ! do ineigh
  end do ! do iatom
  return
end subroutine assemble_ca_2c_dip
