subroutine Dassemble_ca_olsxc_2c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, ratom, imass, neigh_b, neigh_j, neighn, neigh_self, numorb_max, sp_mat, Qin, rho, rho_off, rhoij_off, &
    & s_mat, arhoij_off, arho_off, arhopij_off, arhop_off, rhop_off, rhopij_off, xl, fotxc, fotxc_ca
  use M_fdata, only: num_orb, nssh, lssh, Qneutral
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
  fotxc = 0.0d0
  fotxc_ca = 0.0d0
  do iatom = 1, natoms
   matom = neigh_self(iatom)
   r1(:) = ratom(:,iatom)
   in1 = imass(iatom)
   do ineigh = 1, neighn(iatom) 
    mbeta = neigh_b(ineigh,iatom)
    jatom = neigh_j(ineigh,iatom)
    r2(:) = ratom(:,jatom) + xl(:,mbeta)
    in2 = imass(jatom)
    if (iatom .eq. jatom .and. mbeta .eq. 0) then
    else
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
     interaction = 6
     in3 = in2
     call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcxcx, bcxcpx)
     do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom) - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
       end do
      end do
     end do 
     interaction = 18
     in3 = in2
     do isorp = 1, nssh(in1)
       call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, bcxcx, bcxcpx )
       dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
       do inu = 1, num_orb(in3)
         do imu = 1, num_orb(in1)
           do ix = 1, 3
             fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom) - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)*dxn
           end do
         end do
       end do
     end do
     interaction = 19
     in3 = in2
     do isorp = 1, nssh(in2)
      call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcxcx, bcxcpx )
      dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
      do inu = 1, num_orb(in3)
       do imu = 1, num_orb(in1)
        do ix = 1, 3
         fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom)  - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)*dxn
        end do
       end do
      end do
     end do
     n1 = 0
     do issh = 1, nssh(in1)
      l1 = lssh(issh,in1)
      n1 = n1 + l1 + 1
      n2 = 0
      do jssh = 1, nssh(in2)
       l2 = lssh(jssh,in2)
       n2 = n2 + l2 + 1
       rho_av =  arhoij_off(issh,jssh,ineigh,iatom)
       rhop_av(:) =  arhopij_off(:,issh,jssh,ineigh,iatom)
       call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
       do index1 = -l1, l1
        do index2 = -l2, l2
         imu = n1 + index1
         inu = n2 + index2
         sx = s_mat(imu,inu,ineigh,iatom)
         spx(:) = sp_mat(:,imu,inu,ineigh,iatom)
         rhoin =  rhoij_off(imu,inu,ineigh,iatom)
         rhoinp(:) =  rhopij_off(:,imu,inu,ineigh,iatom)
         mxcb(:,imu,inu) = spx(:)*(rho_av*dmuxc - muxc) + rhop_av(:)*d2muxc*(rho_av*sx - rhoin) - rhoinp(:)*dmuxc
        end do ! do index
       end do ! do index
       n2 = n2 + l2
      end do
      n1 = n1 + l1
     end do
     n1 = 0
     do issh = 1, nssh(in1)
      l1 = lssh(issh,in1)
      n1 = n1 + l1 + 1
      n2 = 0
      do jssh = 1, nssh(in2)
       l2 = lssh(jssh,in2)
       n2 = n2 + l2 + 1
       rho_av =  arho_off(issh,jssh,ineigh,iatom)
       rhop_av =  arhop_off(:,issh,jssh,ineigh,iatom)
       call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
       do index1 = -l1, l1
        do index2 = -l2, l2
         imu = n1 + index1
         inu = n2 + index2 
         sx = s_mat(imu,inu,ineigh,iatom)
         spx(:) = sp_mat(:,imu,inu,ineigh,iatom)
         rhoin =  rho_off(imu,inu,ineigh,iatom)
         rhoinp(:) =  rhop_off(:,imu,inu,ineigh,iatom)
         mxcb(:,imu,inu) = mxcb(:,imu,inu) + spx(:)*(muxc - rho_av*dmuxc)  + rhop_av(:)*d2muxc*( rhoin - rho_av*sx ) + rhoinp(:)*dmuxc
         fotxc_ca(:,ineigh,iatom) = fotxc_ca(:,ineigh,iatom)  - rho(imu,inu,ineigh,iatom)*mxcb(:,imu,inu)
        end do ! do index
       end do ! do index
       n2 = n2 + l2
      end do
      n1 = n1 + l1
     end do
    end if
   end do
  end do
end subroutine Dassemble_ca_olsxc_2c
