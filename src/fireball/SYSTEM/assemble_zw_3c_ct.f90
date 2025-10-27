subroutine assemble_zw_3c_ct (nprocs, iordern, igauss)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, natoms, neigh_comn, neigh_comm, neigh_comj, neigh_comb, neigh_com_ng, neigh_back, imass, ratom, xl, nssh, Qneutral, Qin, orb2shell, s_mat, dip, vxc_ca, g2nu, gvhxc, iqout, Kscf

  implicit none
  integer ialp
  integer iatom
  integer iatomstart
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
  integer kneigh
  integer ineigh1
  integer ineigh2
  integer jbeta
  integer jcount
  integer jcount_sav
  integer jssh
  integer issh1
  integer issh2
  integer matom
  integer mbeta
  integer mneigh
  integer my_proc
  integer natomsp
  integer ix
  integer j
   
  real cost
  real distance_13
  real distance_23
  real dq3
  real dstn_temp
  real dterm
  real dxn
  real rcutoff_ialp
  real rend1
  real rend2
  real sterm
  real stn_temp1
  real stn_temp2
  real x
  real y
  real rcutoff_i
  real rcutoff_j
  real dot_product_dipc_x
  real :: A,B
   
  real, dimension (numorb_max, numorb_max) :: bcca
  real, dimension (numorb_max, numorb_max) :: bccax
  real, dimension (numorb_max, numorb_max) :: emnpl
  real, dimension (3, 3, 3) :: deps
  real, dimension (3, 3) :: eps
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
  real, dimension (3) :: rhat
  real, dimension (3) :: rna
  real, dimension (3) :: rnabc
  real, dimension (3) :: sighat
  
  real stn1
  real stn2
  
  real, dimension (:,:), allocatable :: smG
  real, dimension (:,:,:), allocatable :: spmG
  real, dimension (:, :, :, :), allocatable :: smatG
  real, dimension (:, :, :, :), allocatable :: spmatG
        
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    dq3 = 0.0d0
    do issh = 1, nssh(indna)
      dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
    end do
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh,ialp)
      if (mneigh .ne. 0) then
        iatom = neigh_comj(1,ineigh,ialp)
        ibeta = neigh_comb(1,ineigh,ialp)
        ineigh1 = neigh_com_ng(1,ineigh,ialp)
        r1(:) = ratom(:,iatom) + xl(:,ibeta)
        in1 = imass(iatom)
        jatom = neigh_comj(2,ineigh,ialp)
        jbeta = neigh_comb(2,ineigh,ialp)
        ineigh2 = neigh_com_ng(2,ineigh,ialp)
        r2(:) = ratom(:,jatom) + xl(:,jbeta)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,mneigh)
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
        bcca = 0.0d0
        bccax=0.0d0
        do isorp = 1, nssh(indna)
          dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              issh1=orb2shell(imu,in1)
              issh2=orb2shell(inu,in2)
              A=0.5*s_mat(imu,inu,mneigh,iatom)-dip(imu,inu,mneigh,iatom)/y
              B=0.5*s_mat(imu,inu,mneigh,iatom)+dip(imu,inu,mneigh,iatom)/y  
              bccax(imu,inu) = bccax(imu,inu)+ (A*g2nu(isorp,issh1,ineigh1,ialp)+B*g2nu(isorp,issh2,ineigh2,ialp))*dxn
              if (Kscf .eq. 1 .and. iqout .eq. 6) then 
                gvhxc(imu,inu,isorp,ialp,mneigh,iatom) = gvhxc(imu,inu,isorp,ialp,mneigh,iatom) + A*g2nu(isorp,issh1,ineigh1,ialp)+B*g2nu(isorp,issh2,ineigh2,ialp)
                gvhxc(inu,imu,isorp,ialp,jneigh,jatom) = gvhxc(imu,inu,isorp,ialp,mneigh,iatom)
            end if 
            end do 
          end do 
        end do 
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vxc_ca(imu,inu,mneigh,iatom) = vxc_ca(imu,inu,mneigh,iatom) +  bccax(imu,inu)
            vxc_ca(inu,imu,jneigh,jatom) = vxc_ca(imu,inu,mneigh,iatom)
          end do  
        end do  
      end if
    end do
  end do
  601 format (9(f8.4))
  602 format (9(f8.4))  
  return 
end
