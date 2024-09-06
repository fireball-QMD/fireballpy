subroutine Dassemble_lr_dip ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system
  use M_fdata, only: nssh, Qneutral,num_orb
  implicit none
  integer(c_long) iatom
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) ineigh
  integer(c_long) inu
  integer(c_long) issh
  integer(c_long) jatom
  integer(c_long) ialp
  integer(c_long) inalp
  integer(c_long) ix
  integer(c_long) jneigh
  integer(c_long) ipair
  real(c_double) dist13
  real(c_double) dist23
  real(c_double) dq3
  real(c_double) dterm 
  real(c_double) x
  real(c_double) sterm
  real(c_double), dimension(3)  :: rnabc
  real(c_double), dimension(3)  :: r13
  real(c_double), dimension(3)  :: r23
  real(c_double), dimension(3)  :: r21
  real(c_double), dimension(3)  :: ddterm
  real(c_double), dimension(3)  :: dptermA
  real(c_double), dimension(3)  :: dptermB
  real(c_double), dimension(3)  :: rna
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: spterm
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl
  real(c_double), dimension (3, numorb_max, numorb_max) :: demnplA
  real(c_double), dimension (3, numorb_max, numorb_max) :: demnplB
  real(c_double), dimension (3, numorb_max, numorb_max) :: demnplC

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
end

