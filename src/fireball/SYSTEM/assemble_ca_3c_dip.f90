subroutine assemble_ca_3c_dip ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: eq2
  use M_system, only: smt_elect, natoms, ratom, imass, neigh_max, ewaldsr, neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, &
    & neigh_comm, neigh_comn, neigh_back, neigh_self, numorb_max, Qin, s_mat, vca, dipc, xl, g_h, get_shell_ofatom_issh, kscf, iqout

  use M_fdata, only: nssh, Qneutral,rcutoff, num_orb
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
  integer jatom
  integer jneigh
  integer jbeta
  integer matom
  integer mbeta
  integer mneigh
  integer alpha, beta
  real(dp) aux
  real(dp) cost
  real(dp) distance_13
  real(dp) distance_23
  real(dp) dq3
  real(dp) dstn_temp
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
  real(dp), dimension (numorb_max, numorb_max) :: emnpl
  real(dp), dimension (numorb_max, numorb_max) :: emnpl_noq
  real(dp), dimension (3, 3, 3) :: deps
  real(dp), dimension (3, 3) :: eps
  real(dp), dimension (3) :: r1
  real(dp), dimension (3) :: r2
  real(dp), dimension (3) :: r21
  real(dp), dimension (3) :: rhat
  real(dp), dimension (3) :: rna
  real(dp), dimension (3) :: rnabc
  real(dp), dimension (3) :: sighat
  real(dp) stn1
  real(dp) stn2
  real(dp), dimension (:,:), allocatable :: smG
  real(dp), dimension (:,:,:), allocatable :: spmG
  real(dp), dimension (:, :, :, :), allocatable :: smatG
  real(dp), dimension (:, :, :, :), allocatable :: spmatG
  allocate (smG (numorb_max, numorb_max))
  allocate (spmG (3, numorb_max, numorb_max))
  allocate (smatG (numorb_max, numorb_max, neigh_max, natoms))
  allocate (spmatG (numorb_max, numorb_max, neigh_max, natoms))
  smG = 0.0_dp
  smatG = 0.0_dp
  spmG = 0.0_dp
  spmatG = 0.0_dp
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)  
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      r21(:) = r2(:) - r1(:)
      y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
      if (y .lt. 1.0e-05_dp) then
        sighat(1) = 0.0_dp
        sighat(2) = 0.0_dp
        sighat(3) = 1.0_dp
      else
        sighat(:) = r21(:)/y
      end if
      call epsilon (r2, sighat, eps)
      call deps2cent (r1, r2, eps, deps)
    end do
  end do
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    dq3 = 0.0_dp
    do issh = 1, nssh(indna)
      dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
    end do
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
        rcutoff_i = 0
        do imu = 1, nssh(in1)
          if (rcutoff(imu,in1) .gt. rcutoff_i) rcutoff_i = rcutoff(imu,in1)
        end do
        jatom = neigh_comj(2,ineigh,ialp)
        jbeta = neigh_comb(2,ineigh,ialp)
        r2(:) = ratom(:,jatom) + xl(:,jbeta)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,mneigh)
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
        if (x .lt. 1.0e-05_dp) then
          rhat(1) = 0.0_dp
          rhat(2) = 0.0_dp
          rhat(3) = 0.0_dp
        else
          rhat(:) = rnabc(:)/x
        end if
        cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
        call epsilon (rhat, sighat, eps)
        distance_13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2 + (rna(3) - r1(3))**2)
        distance_23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2 + (rna(3) - r2(3))**2)
        rend1 = rcutoff_i + rcutoff_ialp
        rend2 = rcutoff_j + rcutoff_ialp
        call smoother (distance_13, rend1, smt_elect, stn_temp1, dstn_temp)
        call smoother (distance_23, rend2, smt_elect, stn_temp2, dstn_temp)
        stn1 = stn_temp1*stn_temp2
        stn2 = 1.0_dp - stn1
        if (x .lt. 1.0e-05_dp) then
          emnpl = 0.0_dp
          emnpl_noq = 0.0_dp
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,mneigh,iatom)
              dterm = (dipc(1,imu,inu,mneigh,iatom)*rnabc(1)  + dipc(2,imu,inu,mneigh,iatom)*rnabc(2) + dipc(3,imu,inu,mneigh,iatom)*rnabc(3))
              emnpl(imu,inu) = dq3*sterm/x + dq3*dterm/(x*x*x)           
              emnpl_noq(imu,inu) = sterm/x + dterm/(x*x*x)           
              ewaldsr(imu,inu,mneigh,iatom) = ewaldsr(imu,inu,mneigh,iatom) + emnpl(imu,inu)*eq2
              if (Kscf .eq. 1) then
                do issh = 1, nssh(indna)
                  g_h(get_shell_ofatom_issh(ialp,issh),imu,inu,mneigh,iatom)  =  g_h(get_shell_ofatom_issh(ialp,issh),imu,inu,mneigh,iatom) - emnpl_noq(imu,inu)*eq2
                  g_h(get_shell_ofatom_issh(ialp,issh),inu,imu,jneigh,jatom)  =  g_h(get_shell_ofatom_issh(ialp,issh),imu,inu,mneigh,iatom)
                end do
              end if
              ewaldsr(inu,imu,jneigh,jatom) = ewaldsr(imu,inu,mneigh,iatom)
            end do
          end do
        end if    ! end if (x .lt. 1.0e-05_dp)  
        bcca = 0.0_dp
        do isorp = 1, nssh(indna)
          interaction = 1
          call trescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, bccax)
          dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
              if (Kscf .eq. 1) then
                g_h(get_shell_ofatom_issh(ialp,isorp),imu,inu,mneigh,iatom)  =  g_h(get_shell_ofatom_issh(ialp,isorp),imu,inu,mneigh,iatom) + (stn1*bccax(imu,inu) + stn2*emnpl_noq(imu,inu))*eq2
                g_h(get_shell_ofatom_issh(ialp,isorp),inu,imu,jneigh,jatom)  =  g_h(get_shell_ofatom_issh(ialp,isorp),imu,inu,mneigh,iatom)
              end if
            end do
          end do
        end do ! end do of isorp loop
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vca(imu,inu,mneigh,iatom) = vca(imu,inu,mneigh,iatom) + (stn1*bcca(imu,inu) +  stn2*emnpl(imu,inu))*eq2
            vca(inu,imu,jneigh,jatom) = vca(imu,inu,mneigh,iatom)
          end do
        end do
      end if
    end do
  end do
  deallocate (smG)
  deallocate (spmG)
  deallocate (smatG)
  deallocate (spmatG)
end subroutine assemble_ca_3c_dip
