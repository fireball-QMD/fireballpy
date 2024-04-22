subroutine Dassemble_lr ()
  use M_system
  use M_fdata, only: num_orb,nssh,Qneutral
  use M_constants
  implicit none
  integer iatom
  integer ierror
  integer imu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer inu
  integer issh
  integer jatom
  integer jmu
  integer katom
  integer mbeta
  real(8) distance12
  real(8) dq1
  real(8) dq2
  real(8) dq3
  real(8) dterm
  real(8) sterm
  real(8), dimension (3) :: dewaldlr_i
  real(8), dimension (3) :: dewaldlr_j
  real(8), dimension (3) :: dewaldlr_k
  real(8), dimension (3) :: dpterm
  real(8), dimension (3) :: r1
  real(8), dimension (3) :: r2
  real(8), dimension (3) :: rhat12
  real(8), dimension (3) :: spterm
  real(8), dimension (natoms) :: sub_ewald
  real(8), dimension (3, natoms) :: sub_dewald
  flrew = 0.0d0
  sub_ewald = 0.0d0
  sub_dewald = 0.0d0
  do iatom = 1, natoms
    do jatom = 1, natoms
      in2 = imass(jatom)
      dq2 = 0.0d0
      do issh = 1, nssh(in2)
        dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      sub_ewald(iatom) = sub_ewald(iatom) + dq2*ewald(iatom,jatom)
      sub_dewald(:,iatom) = sub_dewald(:,iatom) + dq2*dewald(:,iatom,jatom)
    end do
  end do
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    dq1 = 0.0d0
    do issh = 1, nssh(in1)
      dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
    end do
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      mbeta = neigh_b(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      dq2 = 0.0d0
      do issh = 1, nssh(in2)
        dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      distance12 = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2 + (r2(3) - r1(3))**2)
      if (distance12 .gt. 1.0d-4) then
        rhat12(:) = (r2(:) - r1(:))/distance12
      end if
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
          if (distance12 .gt. 1.0d-4) then
            dterm = dip(imu,inu,ineigh,iatom)/distance12
            dpterm(:) = dipp(:,imu,inu,ineigh,iatom)/distance12  + dip(imu,inu,ineigh,iatom)*rhat12(:)/distance12**2
            spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)/2.0d0
            dewaldlr_i(:) = (sterm - dterm)*sub_dewald(:,iatom)  + (spterm(:) - dpterm(:))*sub_ewald(iatom)    + (sterm + dterm)*dewald(:,iatom,jatom)*dq1  + (spterm(:) + dpterm(:))*sub_ewald(jatom )
            dewaldlr_j(:) = (sterm - dterm)*dewald(:,jatom,iatom)*dq2  - (spterm(:) - dpterm(:))*sub_ewald(iatom)    + (sterm + dterm)*sub_dewald(:,jatom)  - (spterm(:) + dpterm(:))*sub_ewald(jatom)
          else
            dterm = 0.0d0
            dpterm(:) = 0.0d0
            spterm(:) = 0.0d0
            dewaldlr_i(:) = sterm*sub_dewald(:,iatom)
            dewaldlr_j(:) = sterm*sub_dewald(:,jatom)
          end if
          do jmu = 1, 3
            flrew(jmu,iatom) = flrew(jmu,iatom)  - rho(imu,inu,ineigh,iatom)*dewaldlr_i(jmu)*eq2
            flrew(jmu,jatom) = flrew(jmu,jatom)  - rho(imu,inu,ineigh,iatom)*dewaldlr_j(jmu)*eq2
          end do
          do katom = 1, natoms
            if (katom .ne. iatom .and. katom .ne. jatom) then
              in3 = imass(katom)
              dq3 = 0.0d0
              do issh = 1, nssh(in3)
                dq3 = dq3 + (Qin(issh,katom) - Qneutral(issh,in3))
              end do
              dewaldlr_k(:) = ((sterm-dterm)*dewald(:,katom,iatom) + (sterm+dterm)*dewald(:,katom,jatom))*dq3
              do jmu = 1, 3
                flrew(jmu,katom) = flrew(jmu,katom) - rho(imu,inu,ineigh,iatom)*dewaldlr_k(jmu)*eq2
              end do
            end if
          end do
        end do
      end do
    end do
  end do
end

