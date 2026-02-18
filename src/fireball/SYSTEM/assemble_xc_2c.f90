subroutine assemble_xc_2c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, natoms, neigh_self, ratom, imass, neighn, neigh_b, neigh_j, xl, vxc, rho_off, rhoij_off, s_mat, Kscf, g_xc, Qin
  use M_fdata, only: num_orb, nssh, nsh_max, Qneutral
  implicit none
  integer iatom
  integer iatomstart
  integer imu
  integer in1, in2, in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer jatom
  integer kforce
  integer matom
  integer mbeta
  integer natomsp
  integer issh
  real(double)  y 
  real(double) dxn
  real(double), dimension (numorb_max, numorb_max) :: bcxcx
  real(double), dimension (numorb_max, numorb_max) :: denmx
  real(double), dimension (numorb_max, numorb_max) :: den1x
  real(double), dimension (numorb_max, numorb_max) :: rhomx
  real(double), dimension (3, numorb_max, numorb_max) :: rhompx
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3, 3, 3) :: deps
  real(double), dimension (3) :: r1, r2, r21
  real(double), dimension (3) :: sighat
  real(double), dimension (numorb_max, numorb_max) :: sx
  real(double), dimension (nsh_max) :: dqi
  real(double), dimension (nsh_max) :: dqj


  do iatom = 1, natoms
    matom = neigh_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    dqi = 0.0d0
    do issh = 1, nssh(in1)
      dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
    end do

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
      dqj = 0.0d0
      do issh = 1, nssh(in2)
        dqj(issh) = (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      call epsilon (r2, sighat, eps)
      call deps2cent (r1, r2, eps, deps)
      if (iatom .ne. jatom .or. mbeta .ne. 0) then 
        isorp = 0
        interaction = 6
        in3 = in2
        call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)
          end do
        end do
        interaction = 7
        in3 = in2
        do isorp = 1, nssh(in1)
          call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dqi(isorp)
              if (Kscf .eq. 1) then
                g_xc(imu,inu,isorp,iatom,ineigh,iatom) = g_xc(imu,inu,isorp,iatom,ineigh,iatom) + rhomx(imu,inu)
              end if
            end do
          end do
        end do
        interaction = 8
        in3 = in2
        do isorp = 1, nssh(in2)
          call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dqj(isorp)
              if (Kscf .eq. 1) then
                g_xc(imu,inu,isorp,jatom,ineigh,iatom) = g_xc(imu,inu,isorp,jatom,ineigh,iatom) + rhomx(imu,inu)
              end if
            end do
          end do
        end do
      end if
    end do
  end do
  return
end subroutine assemble_xc_2c
