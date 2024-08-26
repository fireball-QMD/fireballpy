! This routine assembles all of the two-center sVNL (separable pseudopotential) interactions.
subroutine assemble_sVNL ()
  use M_constants, only: wp
  use M_system
  use M_fdata, only: num_orbPP, num_orb
  implicit none
 
  integer iatom
  integer imu
  integer in1
  integer in2
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer jatom
  integer matom
  integer mbeta
 
  real(wp) y
  real(wp), dimension (3, 3) :: eps
  real(wp), dimension (3, 3, 3) :: deps
  real(wp), dimension (3) :: r1
  real(wp), dimension (3) :: r2
  real(wp), dimension (3) :: r21
  real(wp), dimension (3) :: sighat
  real(wp), dimension (numorb_max, numorb_max) :: sVNLx
  real(wp), dimension (3, numorb_max, numorb_max) :: spVNLx
 
  do iatom = 1, natoms 
    matom = nPP_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, nPPn(iatom) 
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
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
    
      isorp = 0
      interaction = 5
      call doscentrosPP (interaction, isorp, y, eps, deps, iforce, in1, in2, sVNLx, spVNLx)
      if (ineigh .ne. matom) then
        do inu = 1, num_orbPP(in2)
          do imu = 1, num_orb(in1)
            sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
            if (iforce .eq. 1) then
              spVNL(:,imu,inu,ineigh,iatom) = spVNLx(:,imu,inu)
            end if
          end do
        end do
      else
        do inu = 1, num_orbPP(in2)
          do imu = 1, num_orb(in1)
            sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
            spVNL(:,imu,inu,ineigh,iatom) = 0.0d0
          end do
        end do
      end if
    end do !ineigh
  end do !iatom
end subroutine assemble_sVNL
