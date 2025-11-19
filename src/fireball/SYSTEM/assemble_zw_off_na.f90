subroutine assemble_zw_off_na ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, natoms, neigh_self, ratom, imass, neighn, neigh_b, neigh_j, xl, vxc, rho_off, rhoij_off, s_mat
  use M_fdata, only: num_orb
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
  real  y 
  real dxn
  real, dimension (numorb_max, numorb_max) :: bcxcx
  real, dimension (numorb_max, numorb_max) :: denmx
  real, dimension (numorb_max, numorb_max) :: den1x
  real, dimension (numorb_max, numorb_max) :: rhomx
  real, dimension (3, numorb_max, numorb_max) :: rhompx
  real, dimension (3, 3) :: eps
  real, dimension (3, 3, 3) :: deps
  real, dimension (3) :: r1, r2, r21
  real, dimension (3) :: sighat
  real, dimension (numorb_max, numorb_max) :: sx
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
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            denmx(imu,inu) = rho_off(imu,inu,ineigh,iatom)   
            den1x(imu,inu) = rhoij_off(imu,inu,ineigh,iatom) 
            sx(imu,inu) = s_mat(imu,inu,ineigh,iatom)
          end do
        end do
        call build_zw_off_na (in1, in2, den1x, denmx, sx, ineigh, iatom, bcxcx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
          end do
        end do
      end if
    end do
  end do
  return
end subroutine assemble_zw_off_na
