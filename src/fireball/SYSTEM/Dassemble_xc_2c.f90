subroutine Dassemble_xc_2c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, ratom, imass, neigh_b, neigh_j, neighn, neigh_self, numorb_max, sp_mat, Qin, rho, rho_off, rhoij_off, &
    & s_mat, arhoij_off, arho_off, arhopij_off, arhop_off, rhop_off, rhopij_off, xl, fotxc_ca
  use M_fdata, only: num_orb, nssh, lssh, Qneutral, nsh_max, TWOCENTER_VXC_0, TWOCENTER_VXC_L, TWOCENTER_VXC_R
  implicit none
  integer iatom
  integer inu
  integer imu
  integer in1
  integer in2
  integer in3
  integer index1
  integer index2
  integer ineigh
  integer interaction
  integer isorp
  integer issh
  integer ix
  integer jatom
  integer jssh
  integer kforce
  integer l1
  integer l2
  integer mbeta
  integer matom
  integer n1
  integer n2
  real(double) y
  real(double) muxc
  real(double) dmuxc
  real(double) d2muxc
  real(double) exc
  real(double) dexc
  real(double) d2exc
  real(double) sx
  real(double) rho_av
  real(double) rhoin
  real(double) dxn
  real(double), dimension (numorb_max, numorb_max) :: bcxcx
  real(double), dimension (3, numorb_max, numorb_max) :: bcxcpx
  real(double), dimension (3, numorb_max, numorb_max) :: mxcb
  real(double), dimension (3) :: rhoinp
  real(double), dimension (3) :: rhop_av
  real(double), dimension (3) :: spx
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3, 3, 3) :: deps
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r21
  real(double), dimension (3) :: sighat
  real(double), dimension (nsh_max) :: dqi
  real(double), dimension (nsh_max) :: dqj
  fotxc_ca = 0.0d0
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
    if (iatom .eq. jatom .and. mbeta .eq. 0) then
    else
     dqj = 0.0d0
     do issh = 1, nssh(in2)
       dqj(issh) = (Qin(issh,jatom) - Qneutral(issh,in2))
     end do
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
     kforce = 1
     isorp = 0
     interaction = TWOCENTER_VXC_0
     in3 = in2
     call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcxcx, bcxcpx)
     do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom) - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
       end do
      end do
     end do 
     interaction = TWOCENTER_VXC_L
     in3 = in2
     do isorp = 1, nssh(in1)
       call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, bcxcx, bcxcpx )
       do inu = 1, num_orb(in3)
         do imu = 1, num_orb(in1)
           do ix = 1, 3
             fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom) - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)*dqi(isorp)
           end do
         end do
       end do
     end do
     interaction = TWOCENTER_VXC_R
     in3 = in2
     do isorp = 1, nssh(in2)
      call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcxcx, bcxcpx )
      do inu = 1, num_orb(in3)
       do imu = 1, num_orb(in1)
        do ix = 1, 3
         fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom)  - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)*dqj(isorp)
        end do
       end do
      end do
     end do
    end if
   end do
  end do
end subroutine Dassemble_xc_2c
