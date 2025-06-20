subroutine assemble_lr () 
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, ewald, dip, neigh_b, neigh_j, neighn, Qin, s_mat, ewaldlr, xl
  use M_fdata, only: nssh,Qneutral,num_orb
  implicit none
  integer iatom
  integer imu
  integer inu
  integer in1
  integer in2
  integer ineigh
  integer issh
  integer jatom
  integer mbeta
  real(double) distance12
  real(double) dq1
  real(double) dterm
  real(double) sterm
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (natoms) :: sub_ewald
  ewaldlr = 0.0d0
  sub_ewald = 0.0d0
  do iatom = 1, natoms
    do jatom = 1, natoms
      in2 = imass(jatom)
      dq1 = 0.0d0
      do issh = 1, nssh(in2)
        dq1 = dq1 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      sub_ewald(iatom) = sub_ewald(iatom) + dq1*ewald(iatom,jatom)
    end do
  end do
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      jatom = neigh_j(ineigh,iatom)
      mbeta = neigh_b(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      distance12 = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2 + (r2(3) - r1(3))**2)
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
          if (distance12 .gt. 1.0d-4) then
            dterm = dip(imu,inu,ineigh,iatom)/distance12
          else
            dterm = 0.0d0
          endif
          ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)  + (sterm - dterm)*sub_ewald(iatom)*eq2  + (sterm + dterm)*sub_ewald(jatom)*eq2
       end do
      end do
    end do
  end do
  return
end subroutine assemble_lr
