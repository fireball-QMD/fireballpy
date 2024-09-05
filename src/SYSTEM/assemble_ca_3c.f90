subroutine assemble_ca_3c ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system
  use M_fdata, only: nssh, Qneutral, rcutoff, lssh, num_orb, isorpmax, nspecies
  implicit none
  integer(c_long) ialp
  integer(c_long) iatom
  integer(c_long) ibeta
  integer(c_long) icount
  integer(c_long) icount_sav
  integer(c_long) ierror
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) indna
  integer(c_long) ineigh
  integer(c_long) interaction
  integer(c_long) inu
  integer(c_long) isorp
  integer(c_long) issh
  integer(c_long) jatom
  integer(c_long) jneigh
  integer(c_long) jbeta
  integer(c_long) jcount
  integer(c_long) jcount_sav
  integer(c_long) jssh
  integer(c_long) matom
  integer(c_long) mbeta
  integer(c_long) mneigh
 
  real(c_double) cost
  real(c_double) distance_13
  real(c_double) distance_23
  real(c_double) dq3
  real(c_double) dstn_temp
  real(c_double) dterm
  real(c_double) dxn
  real(c_double) rcutoff_ialp
  real(c_double) rend1
  real(c_double) rend2
  real(c_double) sterm
  real(c_double) stn_temp1
  real(c_double) stn_temp2
  real(c_double) x
  real(c_double) y
  real(c_double) rcutoff_i
  real(c_double) rcutoff_j
 
  real(c_double), dimension (numorb_max, numorb_max) :: bcca
  real(c_double), dimension (numorb_max, numorb_max) :: bccax
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl_noq
  real(c_double), dimension (3, 3, 3) :: deps
  real(c_double), dimension (3, 3) :: eps
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r21
  real(c_double), dimension (3) :: rhat
  real(c_double), dimension (3) :: rna
  real(c_double), dimension (3) :: rnabc
  real(c_double), dimension (3) :: sighat
  real(c_double), dimension (numorb_max, numorb_max) :: stn1
  real(c_double), dimension (numorb_max, numorb_max) :: stn2

  real(c_double), dimension (:,:), allocatable :: smG
  real(c_double), dimension (:,:,:), allocatable :: spmG
  real(c_double), dimension (:, :, :, :), allocatable :: smatG
  real(c_double), dimension (:, :, :, :), allocatable :: spmatG
         
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
      if (rcutoff(indna,imu) .gt. rcutoff_ialp)              &
       rcutoff_ialp = rcutoff(indna,imu)
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
        distance_13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2  + (rna(3) - r1(3))**2)
        distance_23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2  + (rna(3) - r2(3))**2)

        icount_sav = 0
        do issh = 1, nssh(in1)
          jcount_sav = 0
          do jssh = 1, nssh(in2)
            rend1 = rcutoff_i + rcutoff_ialp
            rend2 = rcutoff_j + rcutoff_ialp
            call smoother (distance_13, rend1, smt_elect, stn_temp1, dstn_temp)
            call smoother (distance_23, rend2, smt_elect, stn_temp2, dstn_temp)
            stn_temp1 = stn_temp1*stn_temp2
            stn_temp2 = 1.0d0 - stn_temp1
            do inu = 1, lssh(issh,in1)*2 + 1
              icount = icount_sav + inu
              do imu = 1, lssh(jssh,in2)*2 + 1
                jcount = jcount_sav + imu
                stn1(icount,jcount) = stn_temp1
                stn2(icount,jcount) = stn_temp2
              end do
            end do
            jcount_sav = jcount
          end do
          icount_sav = icount
        end do
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            sterm = dq3*s_mat(imu,inu,mneigh,iatom)/2.0d0
            dterm = dq3*dip(imu,inu,mneigh,iatom)/y
            emnpl(imu,inu) = (sterm - dterm)/distance_13  + (sterm + dterm)/distance_23
            emnpl_noq(imu,inu) = ((s_mat(imu,inu,mneigh,iatom)/2.0d0)- (dip(imu,inu,mneigh,iatom)/y))/distance_13+  ((s_mat(imu,inu,mneigh,iatom)/2.0d0)+ (dip(imu,inu,mneigh,iatom)/y))/distance_23
            ewaldsr(imu,inu,mneigh,iatom) =  ewaldsr(imu,inu,mneigh,iatom) + emnpl(imu,inu)*eq2
            ewaldsr(inu,imu,jneigh,jatom)=ewaldsr(imu,inu,mneigh,iatom)
          end do
        end do
        !  CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
        bcca = 0.0d0
        do isorp = 1, nssh(indna)
          interaction = 1
          call trescentros (interaction, isorp, isorpmax, in1, in2,    indna, x, y, cost, eps, bccax)
          dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
            end do
          end do
        end do ! end do of isorp loop
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vca(imu,inu,mneigh,iatom) = vca(imu,inu,mneigh,iatom) +    (stn1(imu,inu)*bcca(imu,inu) + stn2(imu,inu)*emnpl(imu,inu))*eq2
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
  return
end
