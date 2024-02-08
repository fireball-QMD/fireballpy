! This routine assembles all of the two-center and degenerate two-center interactions.
subroutine assemble_2c ()
  use M_system
  use M_fdata, only: num_orb
  integer iatom
  integer iatomstart
  integer ierror
  integer imu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer jatom
  integer kforce
  integer matom
  integer mbeta
  integer mneigh_self
  integer my_proc
  integer natomsp
  integer ix
  integer iy
  integer iz
 
  real y
  real, dimension (numorb_max, numorb_max) :: bcna
  real, dimension (3, numorb_max, numorb_max) :: bcnapx
  real, dimension (numorb_max, numorb_max) :: bcnax
  real, dimension (3, 3, 3) :: deps
  real, dimension (3, 3) :: eps
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
  real, dimension (3) :: sighat
  real, dimension (numorb_max, numorb_max) :: sx
  real, dimension (3, numorb_max, numorb_max) :: spx
  real, dimension (numorb_max, numorb_max) :: tx
  real, dimension (3, numorb_max, numorb_max) :: tpx
  real, dimension (numorb_max, numorb_max) :: dipx
  real, dimension (3, numorb_max, numorb_max) :: dippx

  vna = 0.0d0
  s_mat = 0.0d0
  t_mat = 0.0d0
  sp_mat = 0.0d0
  tp_mat = 0.0d0
  dipcm = 0.0d0
  dipc = 0.0d0
  dippcm = 0.0d0
  dippc = 0.0d0

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
 
      ! CALL DOSCENTROS AND GET S AND T
      ! ****************************************************************************
      isorp = 0
      interaction = 1
      in3 = in2
      call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, sx, spx)
 
      isorp = 0
      interaction = 13
      in3 = in2
      call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, tx, tpx)
 
      ! Write s and t to appropriate arrays
      do inu = 1, num_orb(in2)
        do imu = 1, num_orb(in1)
          s_mat(imu,inu,ineigh,iatom) = sx(imu,inu)
          t_mat(imu,inu,ineigh,iatom) = tx(imu,inu)
          if (iforce .eq. 1) then
            sp_mat(:,imu,inu,ineigh,iatom) = spx(:,imu,inu)
            tp_mat(:,imu,inu,ineigh,iatom) = tpx(:,imu,inu)
          end if
        end do
      end do

      if (ineigh .eq. matom) then
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in2)
            s_mat(imu,inu,ineigh,iatom) = 0.0d0
          end do
          s_mat(inu,inu,ineigh,iatom) = 1.0d0
        end do
      end if

      isorp = 0
      kforce = 1           ! don't calculate forces here
      interaction = 4
      in3 = in1
      call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcnax, bcnapx)
 
      do inu = 1, num_orb(in3)
        do imu = 1, num_orb(in1)
          vna(imu,inu,matom,iatom) =  vna(imu,inu,matom,iatom) + bcnax(imu,inu)*eq2
        end do
      end do


      ! CALL DOSCENTROS AND GET VNA FOR ONTOP CASE
      ! ****************************************************************************
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        ! Do nothing here - special case. Interaction already calculated in atm case.
      else
        isorp = 0
        interaction = 2
        in3 = in2
        call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, bcnax, bcnapx)
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            bcna(imu,inu) = bcnax(imu,inu)
          end do
        end do
        ! For the vna_ontopr case, the potential is in the second atom (jatom): Neutral atom piece
        isorp = 0
        interaction = 3
        in3 = in2
        call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, bcnax, bcnapx)
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            bcna(imu,inu) = bcna(imu,inu) + bcnax(imu,inu)
          end do
        end do
        ! Now put into vna.
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            vna(imu,inu,ineigh,iatom) =  vna(imu,inu,ineigh,iatom) + bcna(imu,inu)*eq2
          end do
        end do
      end if ! End if for r1 .ne. r2 case


      ! JIMM: we read here the Z,Y,X dipole matrix elements and derivatives for the dipole long-range theory
      if (idipole .eq. 1) then
        ! CALL DOSCENTROS AND GET DIP Z
        isorp = 0
        interaction = 9
        in3 = in2
        call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, dipx, dippx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            dip(imu,inu,ineigh,iatom) = dipx(imu,inu)
            dipcm(3,imu,inu) = dipx(imu,inu)
            if (iforce .eq. 1) then 
              dippcm(:,3,imu,inu) = dippx(:,imu,inu)
            end if ! end if iforce = 1
          end do
        end do


        ! CALL DOSCENTROS AND GET DIP Y
        isorp = 0
        interaction = 10
        in3 = in2
        call doscentrosDipY (interaction, isorp, in1, in2, in3, y, eps, deps, dipx, dippx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            dipcm(2,imu,inu) = dipx(imu,inu)
            if (iforce .eq. 1) dippcm(:,2,imu,inu) = dippx(:,imu,inu)
          end do
        end do

        ! CALL DOSCENTROS AND GET DIP X
        isorp = 0
        interaction = 11
        in3 = in2
        call doscentrosDipX (interaction, isorp, in1, in2, in3, y, eps, deps, dipx, dippx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            dipcm(1,imu,inu) = dipx(imu,inu)
            if (iforce .eq. 1) dippcm(:,1,imu,inu)  = dippx(:,imu,inu)
          end do
        end do
        do ix = 1, 3
          do iy = 1, 3
            dipc(ix,:,:,ineigh,iatom) = dipc(ix,:,:,ineigh,iatom) + eps(ix,iy) *dipcm(iy,:,:)
          enddo
        enddo
        do ix = 1,3
          do iy = 1,3
            do iz = 1,3
              dippc(ix,iy,:,:,ineigh,iatom) = dippc(ix,iy,:,:,ineigh,iatom)+deps(ix,iy,iz)*dipcm(iz,:,:) + eps(iy,iz)*dippcm(ix,iz,:,:)
            enddo
          enddo
        enddo
      end if ! idipole = 1
    end do
  end do
return
end subroutine assemble_2c
