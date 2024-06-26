subroutine assemble_ca_3c_dip ()
  use M_system
  use M_constants
  use M_fdata, only: nssh, Qneutral,rcutoff, num_orb, isorpmax, nspecies
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer icount
  integer icount_sav
  integer ierror
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
  integer jcount
  integer jcount_sav
  integer jssh
  integer matom
  integer mbeta
  integer mneigh
  integer ix
  integer j
  real*8 cost
  real*8 distance_13
  real*8 distance_23
  real*8 dq3
  real*8 dstn_temp
  real*8 dterm
  real*8 dxn
  real*8 rcutoff_ialp
  real*8 rend1
  real*8 rend2
  real*8 sterm
  real*8 stn_temp1
  real*8 stn_temp2
  real*8 x
  real*8 y
  real*8 rcutoff_i
  real*8 rcutoff_j
  real*8 dot_product_dipc_x
  real*8, dimension (numorb_max, numorb_max) :: bcca
  real*8, dimension (numorb_max, numorb_max) :: bccax
  real*8, dimension (numorb_max, numorb_max) :: emnpl
  real*8, dimension (numorb_max, numorb_max) :: emnpl_noq
  real*8, dimension (3, 3, 3) :: deps
  real*8, dimension (3, 3) :: eps
  real*8, dimension (3) :: r1
  real*8, dimension (3) :: r2
  real*8, dimension (3) :: r21
  real*8, dimension (3) :: rhat
  real*8, dimension (3) :: rna
  real*8, dimension (3) :: rnabc
  real*8, dimension (3) :: sighat
  real*8 stn1
  real*8 stn2
  real*8, dimension (:,:), allocatable :: smG
  real*8, dimension (:,:,:), allocatable :: spmG
  real*8, dimension (:, :, :, :), allocatable :: smatG
  real*8, dimension (:, :, :, :), allocatable :: spmatG
  allocate (smG (numorb_max, numorb_max))
  allocate (spmG (3, numorb_max, numorb_max))
  allocate (smatG (numorb_max, numorb_max, neigh_max, natoms))
  allocate (spmatG (numorb_max, numorb_max, neigh_max, natoms))
  smG = 0.0d0
  smatG = 0.0d0
  spmG = 0.0d0
  spmatG = 0.0d0
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
      if (y .lt. 1.0d-05) then
        sighat(1) = 0.0d0
        sighat(2) = 0.0d0
        sighat(3) = 1.0d0
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
    dq3 = 0.0d0
    do issh = 1, nssh(indna)
      dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
    end do
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
        rcutoff_i = 0
        do imu = 1, nssh(in1)
          if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
        end do
        jatom = neigh_comj(2,ineigh,ialp)
        jbeta = neigh_comb(2,ineigh,ialp)
        r2(:) = ratom(:,jatom) + xl(:,jbeta)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,mneigh)
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
        if (x .lt. 1.0d-05) then
          rhat(1) = 0.0d0
          rhat(2) = 0.0d0
          rhat(3) = 0.0d0
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
        stn2 = 1.0d0 - stn1
        if (x .lt. 1.0d-05) then
          emnpl = 0.0d0
          emnpl_noq = 0.0d0
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,mneigh,iatom)
              dterm = (dipc(1,imu,inu,mneigh,iatom)*rnabc(1)  + dipc(2,imu,inu,mneigh,iatom)*rnabc(2) + dipc(3,imu,inu,mneigh,iatom)*rnabc(3))
              emnpl(imu,inu) = dq3*sterm/x + dq3*dterm/(x*x*x)           
              emnpl_noq(imu,inu) = sterm/x + dterm/(x*x*x)           
              ewaldsr(imu,inu,mneigh,iatom) = ewaldsr(imu,inu,mneigh,iatom) + emnpl(imu,inu)*eq2
              ewaldsr(inu,imu,jneigh,jatom) = ewaldsr(imu,inu,mneigh,iatom)
            end do
          end do
        end if    ! end if (x .lt. 1.0d-05)  
        bcca = 0.0d0
        do isorp = 1, nssh(indna)
          interaction = 1
          call trescentros (interaction, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, bccax)
          dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
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
end

