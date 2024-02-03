subroutine assemble_lr () 
  use M_system
  use M_constants
  use M_fdata, only: nssh,Qneutral,num_orb
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
  integer natomsp
  integer katom
  integer in3
  real distance12
  real dq1
  real dterm
  real sterm
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (natoms) :: sub_ewald
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
          if (Kscf .eq. 1 .and. iqout .eq. 6) then
            do katom = 1,natoms
              in3 = imass(katom)
              do issh = 1, nssh(in3)
                gvhxc(imu,inu,issh,katom,ineigh,iatom) =  gvhxc(imu,inu,issh,katom,ineigh,iatom) +  (sterm - dterm)*ewald(iatom,katom)*eq2 +  (sterm + dterm)*ewald(jatom,katom)*eq2
              end do ! end do issh
            end do    ! end do katom
          end if !  end if Kscf .eq. 1 .and. iqout .eq. 6
        end do
      end do
    end do
  end do
  return
  end

