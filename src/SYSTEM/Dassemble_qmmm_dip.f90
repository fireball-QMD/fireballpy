subroutine Dassemble_qmmm_dip () 
  use M_system
  use M_fdata, only  : nssh, Qneutral, num_orb
  use M_constants, only : eq2
  implicit none
  integer iatom
  integer imu
  integer inu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer issh
  integer jatom
  integer katom
  integer mbeta
  integer ix
  real dist13
  real dist23
  real distance12
  real dij
  real dq1
  real dq2
  real dq3
  real dq4
  real dterm
  real x
  real sterm
  real, dimension (3) :: rna
  real, dimension (3) :: rnabc
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
  real, dimension (3) :: vij
  real, dimension(3)  :: ddterm
  real, dimension(3)  :: dptermA
  real, dimension(3)  :: dptermB
  real, dimension (3) :: spterm
  real, dimension (numorb_max, numorb_max) :: emnpl
  real, dimension (3, numorb_max, numorb_max) :: demnplA
  real, dimension (3, numorb_max, numorb_max) :: demnplB
  real, dimension (3, numorb_max, numorb_max) :: demnplC
  flrew_qmmm = 0.0d0
  qmmm_dxyzcl = 0.0d0
  do iatom=1, natoms
   r1(:) = ratom(:,iatom)
   in1 = imass(iatom)
   do ineigh = 1, neighn(iatom)
    mbeta = neigh_b(ineigh,iatom)
    jatom = neigh_j(ineigh,iatom)
    r2(:) = ratom(:,jatom) + xl(:,mbeta)
    in2 = imass(jatom)
    do katom = 1, qmmm_qm_mm_pairs
     rna(1) = qmmm_qm_xcrd(1,katom)
     rna(2) = qmmm_qm_xcrd(2,katom)
     rna(3) = qmmm_qm_xcrd(3,katom)
     dq3 = - qmmm_qm_xcrd(4,katom) ! charge in amber have opposite sign
     r21(:) = r2(:) - r1(:)
     rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
     x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       sterm = s_mat(imu,inu,ineigh,iatom)
       dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1)   &
    &   + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)   &
    &   + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))
       ddterm(:) = (dippc(:,1,imu,inu,ineigh,iatom)*rnabc(1)   &
         &  + dippc(:,2,imu,inu,ineigh,iatom)*rnabc(2)   &
         &  + dippc(:,3,imu,inu,ineigh,iatom)*rnabc(3))
       spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)
       dptermA(:)=   dipc(:,imu,inu,ineigh,iatom)/(x*x*x)   &
         &   - 3*dterm*rnabc(:)/(x*x*x*x*x)
       dptermB(:)= - 0.50*dipc(:,imu,inu,ineigh,iatom)/(x*x*x)  &
         &   + 0.50*3*dterm*rnabc(:)/(x*x*x*x*x)    &
         &   + ddterm(:)/(x*x*x)
       demnplA(:,imu,inu) = - dq3*sterm*rnabc(:)/(x*x*x)   &
         &      + dq3*dptermA(:)
       demnplB(:,imu,inu) =  dq3*0.50*sterm*rnabc(:)/(x*x*x)  &
         &     + dq3*spterm(:)/(x)    &
         &     + dq3*dptermB(:)
       demnplC(:,imu,inu) = - demnplA(:,imu,inu) - demnplB(:,imu,inu)
       do ix = 1, 3
        qmmm_dxyzcl(ix,katom)  = qmmm_dxyzcl(ix,katom) &
      &              +rho(imu,inu,ineigh,iatom)*demnplA(ix,imu,inu)*eq2*23.061d0
        flrew(ix,iatom) =  flrew(ix,iatom) &
      &      - rho(imu,inu,ineigh,iatom)*demnplB(ix,imu,inu)*eq2
        flrew(ix,jatom) =  flrew(ix,jatom) &
      &      - rho(imu,inu,ineigh,iatom)*demnplC(ix,imu,inu)*eq2
       end do ! do ix
      end do !end do imu = 1, num_orb(in1)
     end do  ! end do inu = 1, num_orb(in2)
    end do    ! end do katom: mm atom
   end do  ! end do ineigh = 1, neighn(iatom)
  end do   ! end do iatom = 1, natoms
  do iatom=1, natoms
   in3 = imass(iatom)
   dq4 = 0.0d0
   do issh = 1, nssh(in3)
    dq4 = dq4  + Qneutral(issh,in3)
   end do
   do katom = 1, qmmm_qm_mm_pairs
    dij = sqrt ((qmmm_qm_xcrd(1,katom)-ratom(1,iatom))*(qmmm_qm_xcrd(1,katom)-ratom(1,iatom)) + &
          (qmmm_qm_xcrd(2,katom)-ratom(2,iatom))*(qmmm_qm_xcrd(2,katom)-ratom(2,iatom)) + & 
          (qmmm_qm_xcrd(3,katom)-ratom(3,iatom))*(qmmm_qm_xcrd(3,katom)-ratom(3,iatom)) )
    vij(:) = qmmm_qm_xcrd(:,katom)-ratom(:,iatom)
    do ix = 1, 3
     qmmm_dxyzcl(ix,katom) = qmmm_dxyzcl(ix,katom) - dq4*qmmm_qm_xcrd(4,katom)*(vij(ix) / dij**3)*eq2*23.061d0
     flrew_qmmm(ix,iatom) = flrew_qmmm(ix,iatom) - dq4*qmmm_qm_xcrd(4,katom)*(vij(ix) / dij**3)*eq2
    enddo
   end do
  end do
  return
  end

