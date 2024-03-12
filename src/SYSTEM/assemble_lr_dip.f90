subroutine assemble_lr_dip () 
  use M_system
  use M_constants
  use M_fdata, only: nssh, Qneutral, num_orb
  implicit none
  integer iatom
  integer ierror
  integer imu
  integer inu
  integer in1
  integer in2
  integer ineigh
  integer issh
  integer jatom
  integer mbeta
  integer my_proc
  integer ialp
  integer inalp
  integer jneigh
  integer ipair
  real dist13
  real dist23
  real dq3
  real dterm
  real sterm
  real x
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: rna
  real, dimension (3) :: r13
  real, dimension (3) :: r23
  real, dimension (3) :: r21
  real, dimension (3) :: rnabc
  real, dimension (natoms) :: sub_ewald
  real, dimension (numorb_max, numorb_max) :: emnpl
  real, dimension (numorb_max, numorb_max) :: emnpl_noq
  ewaldlr = 0.0d0
  do ipair = 1,tot_pairs
    iatom = neigh_pair_a1(ipair)
    jatom = neigh_pair_a2(ipair)
    ineigh = neigh_pair_n1(ipair)
    jneigh = neigh_pair_n2(ipair)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    r2(:) = ratom(:,jatom)
    in2 = imass(jatom)
    do ialp = 1, natoms    !the ialp is the "distant" atom.
      rna(:) = ratom(:,ialp)
      inalp = imass(ialp)
      dq3 = 0.0d0
      do issh = 1, nssh(inalp)
        dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,inalp))
      end do ! end do issh = 1m nssh(inalp)
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
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,ineigh,iatom)
              dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1) + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)  + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))
              emnpl(imu,inu) = dq3*sterm/x + dq3*dterm/(x*x*x)
              if (Kscf .eq. 1 .and. iqout .eq. 6) then
                emnpl_noq(imu,inu) = sterm/x + dterm/(x*x*x)
              end if 
              ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)  + emnpl(imu,inu)*eq2
              if (Kscf .eq. 1 .and. iqout .eq. 6) then
                do issh = 1, nssh(inalp)
                  gvhxc(imu,inu,issh,ialp,ineigh,iatom) =  gvhxc(imu,inu,issh,ialp,ineigh,iatom) + emnpl_noq(imu,inu)*eq2
                  gvhxc(inu,imu,issh,ialp,jneigh,jatom) =  gvhxc(imu,inu,issh,ialp,ineigh,iatom)  
                end do ! end do issh
              end if ! end if Kscf .eq. 1 iqout .eq. 6
              ewaldlr(inu,imu,jneigh,jatom)  = ewaldlr(imu,inu,ineigh,iatom)
            end do !end do imu = 1, num_orb(in1)
          end do  ! end do inu = 1, num_orb(in2)
        end if    ! end if (x .lt. 1.0d-05)
      end if    ! (dist13 .lt. 1.0d-5) .or. (dist23 .lt. 1.0d-5)
    end do     ! end do ialp = 1, natoms
  end do !end do ipair
  do iatom = 1,natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    jatom = iatom
    r2(:) = ratom(:,jatom)
    in2 = imass(jatom)
    ineigh = neigh_self(iatom)
    do ialp = 1, natoms    !the ialp is the "distant" atom.
      rna(:) = ratom(:,ialp)
      inalp = imass(ialp)
      dq3 = 0.0d0
      do issh = 1, nssh(inalp)
        dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,inalp))
      end do ! end do issh = 1m nssh(inalp)
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
        else
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              sterm = s_mat(imu,inu,ineigh,iatom)
              dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1) + dipc(2,imu,inu,ineigh,iatom)*rnabc(2) + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))
              emnpl(imu,inu) = dq3*sterm/x + dq3*dterm/(x*x*x)
              emnpl_noq(imu,inu) = sterm/x + dterm/(x*x*x)
              ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)  + emnpl(imu,inu)*eq2
              if (Kscf .eq. 1 .and. iqout .eq. 6) then
                do issh = 1, nssh(inalp)
                  gvhxc(imu,inu,issh,ialp,ineigh,iatom) =  gvhxc(imu,inu,issh,ialp,ineigh,iatom)  + emnpl_noq(imu,inu)*eq2
                end do ! end do issh
              end if ! end if Kscf .eq. 1 iqout .eq. 6
            end do !end do imu = 1, num_orb(in1)
          end do  ! end do inu = 1, num_orb(in2)
       end if    ! end if (x .lt. 1.0d-05)
     end if    ! (dist13 .lt. 1.0d-5) .or. (dist23 .lt. 1.0d-5)
   end do     ! end do ialp = 1, natoms
 end do      ! end do iatom = 1,natoms
 return
end

