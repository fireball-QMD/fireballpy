subroutine Dassemble_lr_dip ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, neigh_self, neigh_pair_a1, neigh_pair_a2, neigh_pair_n1, neigh_pair_n2, &
    & numorb_max, tot_pairs, sp_mat, dippc, Qin, rho, s_mat, dipc, flrew
  use M_fdata, only: nssh, Qneutral,num_orb
  implicit none
  integer iatom
  integer imu
  integer in1
  integer in2
  integer ineigh
  integer inu
  integer issh
  integer jatom
  integer ialp
  integer inalp
  integer ix
  integer jneigh
  integer ipair
  real(double) dist13
  real(double) dist23
  real(double) dq3
  real(double) dterm 
  real(double) x
  real(double) sterm
  real(double), dimension(3)  :: rnabc
  real(double), dimension(3)  :: r13
  real(double), dimension(3)  :: r23
  real(double), dimension(3)  :: r21
  real(double), dimension(3)  :: ddterm
  real(double), dimension(3)  :: dptermA
  real(double), dimension(3)  :: dptermB
  real(double), dimension(3)  :: rna
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: spterm
  real(double), dimension (numorb_max, numorb_max) :: emnpl
  real(double), dimension (3, numorb_max, numorb_max) :: demnplA
  real(double), dimension (3, numorb_max, numorb_max) :: demnplB
  real(double), dimension (3, numorb_max, numorb_max) :: demnplC

  flrew = 0.0d0
  do ipair = 1,tot_pairs
    iatom = neigh_pair_a1(ipair)
    jatom = neigh_pair_a2(ipair)
    ineigh = neigh_pair_n1(ipair)
    jneigh = neigh_pair_n2(ipair)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    r2(:) = ratom(:,jatom)
    in2 = imass(jatom)
    do ialp = 1, natoms   !the ialp is the "distant" atom.
      rna(:) = ratom(:,ialp)
      inalp = imass(ialp)
      dq3 = 0.0d0
      do issh = 1, nssh(inalp)
        dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,inalp))
      end do ! end do issh = 1m nssh(inalp)
      in2 = imass(jatom)
      r13=rna(:)-r1(:)
      r23=rna(:)-r2(:)
      dist13=sqrt(r13(1)*r13(1) + r13(2)*r13(2) + r13(3)*r13(3))
      dist23=sqrt(r23(1)*r23(1) + r23(2)*r23(2) + r23(3)*r23(3))
      if ((dist13 .lt. 1.0d-5) .or. (dist23 .lt. 1.0d-5)) then
      else
        r21(:) = r2(:) - r1(:)
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        if (x .lt. 1.0d-05) then
          emnpl = 0.0d0
          demnplA = 0.0d0
          demnplB = 0.0d0
          demnplC = 0.0d0
          else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,ineigh,iatom)
              dterm     = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1)      + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)    + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))
              ddterm(:) = (dippc(:,1,imu,inu,ineigh,iatom)*rnabc(1)   + dippc(:,2,imu,inu,ineigh,iatom)*rnabc(2) + dippc(:,3,imu,inu,ineigh,iatom)*rnabc(3))
              spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)
              dptermA(:)=   dipc(:,imu,inu,ineigh,iatom)/(x*x*x)  - 3*dterm*rnabc(:)/(x*x*x*x*x)
              dptermB(:)= - 0.50*dipc(:,imu,inu,ineigh,iatom)/(x*x*x)  + 0.50*3*dterm*rnabc(:)/(x*x*x*x*x) + ddterm(:)/(x*x*x)
              demnplA(:,imu,inu) = - dq3*sterm*rnabc(:)/(x*x*x) + dq3*dptermA(:)
              demnplB(:,imu,inu) =  dq3*0.50*sterm*rnabc(:)/(x*x*x)  + dq3*spterm(:)/(x)  + dq3*dptermB(:)
              demnplC(:,imu,inu) = - demnplA(:,imu,inu) -demnplB(:,imu,inu)
              do ix = 1, 3
                flrew(ix,ialp)  = flrew(ix,ialp) -2*rho(imu,inu,ineigh,iatom)*demnplA(ix,imu,inu)*eq2
                flrew(ix,iatom) = flrew(ix,iatom) -2*rho(imu,inu,ineigh,iatom)*demnplB(ix,imu,inu)*eq2
                flrew(ix,jatom) = flrew(ix,jatom) -2*rho(imu,inu,ineigh,iatom)*demnplC(ix,imu,inu)*eq2
              end do 
            end do  !imu
          end do   !inu
        end if   ! (x .lt. 1.0d-05)
      end if  ! end if if ((dist13 .lt. 1.0d-5) .or. (dist23 .lt.1.0d-5)) 
    end do     ! end do ialp = 1, natoms
  end do !end do ipair 

  do iatom = 1,natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    ineigh=neigh_self(iatom)
    jatom = iatom
    r2(:) = ratom(:,jatom)
    in2 = imass(jatom)
    do ialp = 1, natoms   !the ialp is the "distant" atom.
      rna(:) = ratom(:,ialp)
      inalp = imass(ialp)
      dq3 = 0.0d0
      do issh = 1, nssh(inalp)
        dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,inalp))
      end do ! end do issh = 1m nssh(inalp)
      in2 = imass(jatom)
      r13=rna(:)-r1(:)
      r23=rna(:)-r2(:)
      dist13=sqrt(r13(1)*r13(1) + r13(2)*r13(2) + r13(3)*r13(3))
      dist23=sqrt(r23(1)*r23(1) + r23(2)*r23(2) + r23(3)*r23(3))
      if ((dist13 .lt. 1.0d-5) .or. (dist23 .lt. 1.0d-5)) then
      else
        r21(:) = r2(:) - r1(:)
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        if (x .lt. 1.0d-05) then
          emnpl = 0.0d0
          demnplA = 0.0d0
          demnplB = 0.0d0
          demnplC = 0.0d0
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,ineigh,iatom)
              dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1)  + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)  + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))
              ddterm(:) = (dippc(:,1,imu,inu,ineigh,iatom)*rnabc(1)  + dippc(:,2,imu,inu,ineigh,iatom)*rnabc(2)  + dippc(:,3,imu,inu,ineigh,iatom)*rnabc(3))
              spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)
              dptermA(:)=   dipc(:,imu,inu,ineigh,iatom)/(x*x*x)  - 3*dterm*rnabc(:)/(x*x*x*x*x)
              dptermB(:)= - 0.50*dipc(:,imu,inu,ineigh,iatom)/(x*x*x) + 0.50*3*dterm*rnabc(:)/(x*x*x*x*x) + ddterm(:)/(x*x*x)
              demnplA(:,imu,inu) = - dq3*sterm*rnabc(:)/(x*x*x)  + dq3*dptermA(:)
              demnplB(:,imu,inu) =  dq3*0.50*sterm*rnabc(:)/(x*x*x)  + dq3*spterm(:)/(x)  + dq3*dptermB(:)
              demnplC(:,imu,inu) = - demnplA(:,imu,inu)-demnplB(:,imu,inu)
              do ix = 1, 3
                flrew(ix,ialp)  = flrew(ix,ialp)  - rho(imu,inu,ineigh,iatom)*demnplA(ix,imu,inu)*eq2
                flrew(ix,iatom) = flrew(ix,iatom) - rho(imu,inu,ineigh,iatom)*demnplB(ix,imu,inu)*eq2
                flrew(ix,jatom) = flrew(ix,jatom) - rho(imu,inu,ineigh,iatom)*demnplC(ix,imu,inu)*eq2
              end do 
            end do  
          end do   
        end if   !(x .lt. 1.0d-05)
      end if  !end if if ((dist13 .lt. 1.0d-5) .or. (dist23 .lt.1.0d-5)) 
    end do !ialp 
  end do !iatom 
end subroutine Dassemble_lr_dip
