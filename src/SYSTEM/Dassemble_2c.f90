subroutine Dassemble_2c ()
  use M_constants, only: wp, eq2
  use M_system        
  use M_fdata, only: num_orb
  implicit none
  integer iatom
  integer ierror
  integer imu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer ix
  integer jatom
  integer kforce
  integer matom
  integer mbeta
  real(wp) muxc
  real(wp) sumS
  real(wp) sumT
  real(wp) y
  real(wp), dimension (numorb_max, numorb_max) :: bcnax
  real(wp), dimension (3, numorb_max, numorb_max) :: bcnapx
  real(wp), dimension (3, 3) :: eps
  real(wp), dimension (3, 3, 3) :: deps
  real(wp), dimension (3) :: r1
  real(wp), dimension (3) :: r2
  real(wp), dimension (3) :: r21
  real(wp), dimension (3) :: sighat
  fana = 0.0d0
  fotna = 0.0d0
  ft = 0.0d0
  fro = 0.0d0
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
      do ix = 1, 3
        sumT = 0.0d0
        sumS = 0.0d0
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            sumT = sumT  + rho(imu,inu,ineigh,iatom)*tp_mat(ix,imu,inu,ineigh,iatom)
            sumS = sumS  + cape(imu,inu,ineigh,iatom)*sp_mat(ix,imu,inu,ineigh,iatom)
          end do
         end do
      ft(ix,iatom) = ft(ix,iatom) + (-1.0d0)*sumT
      fro(ix,iatom) = fro(ix,iatom) +  sumS
      ft(ix,jatom) = ft(ix,jatom) - (-1.0d0)*sumT
      fro(ix,jatom) = fro(ix,jatom) -      sumS
      end do ! do ix
      isorp = 0
      kforce = 1  
      interaction = 4
      in3 = in1
      call doscentros (interaction, isorp, kforce, in1, in2, in3, y,eps, deps, bcnax, bcnapx)
      do inu = 1, num_orb(in3)
       do imu = 1, num_orb(in1)
         do ix = 1, 3
          fana(ix,ineigh,iatom) = fana(ix,ineigh,iatom) - rho(imu,inu,matom,iatom)*bcnapx(ix,imu,inu)*eq2
         end do
       end do
      end do
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
      else
       isorp = 0
       interaction = 2
       in3 = in2
       call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, bcnax, bcnapx)
       do inu = 1, num_orb(in3)
         do imu = 1, num_orb(in1)
          do ix = 1, 3
            fotna(ix,ineigh,iatom) = fotna(ix,ineigh,iatom) - rho(imu,inu,ineigh,iatom)*bcnapx(ix,imu,inu)*eq2
          end do
         end do
       end do
      end if
    end do
  end do
end subroutine Dassemble_2c

