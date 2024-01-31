! This routine assembles all of the two-center and degenerate two-center interactions.
subroutine assemble_2c (iforce)
  use M_system
  integer, intent (in) :: iforce
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
 
      ! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
      ! This gives us the distance dbc (or y value in the 2D grid).
      r21(:) = r2(:) - r1(:)
      y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
      ! Find the unit vector in sigma direction.
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
      ! The derivatives are tpx and spx, where p means derivative and x means crytal
      ! coordinates. The derivatives are:  d/d(b1) <phi(r-b1) ! O ! phi(r-b2)>
      ! where O = 1 (for overlap), or T.  The derivative is a vector in crystal
      ! coordinates and is stored in tp(ix,imu,inu,iatom,ineigh). The subroutine
      ! returns the derivative for just that one value of iatom and ineigh, and the
      ! result is returned in the arguement list, tpx(3,4,4) or spx(3,4,4).  The
      ! derivatives are computed only if iforce = 1 !
      ! For these interactions, there are no subtypes and isorp = 0
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

      ! Sometimes while running diagnostics, we may want to turn off the overlap.
      ! We certainly do not want all of the overlap interactions to be zero, but
      ! instead we would want the unity matrix.
      if (ineigh .eq. matom) then
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in2)
            s_mat(imu,inu,ineigh,iatom) = 0.0d0
          end do
          s_mat(inu,inu,ineigh,iatom) = 1.0d0
        end do
      end if

      ! CALL DOSCENTROS AND GET VNA FOR ATOM CASE
      ! The vna 2 centers are: ontop (L), ontop (R), and atm.
      ! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
      ! If r1 = r2, then this is a case where the two wavefunctions are at the
      ! same site, but the potential vna is at a different site (atm case).
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
      ! if r1 .ne. r2, then this is a case where the potential is located
      ! at one of the sites of a wavefunction (ontop case).
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        ! Do nothing here - special case. Interaction already calculated in atm case.
      else
        ! For the vna_ontopl case, the potential is in the first atom (iatom): Neutral atom piece
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
        call doscentrosDipY (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, dipx, dippx)
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
        call doscentrosDipX (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, dipx, dippx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            dipcm(1,imu,inu) = dipx(imu,inu)
            if (iforce .eq. 1) dippcm(:,1,imu,inu)  = dippx(:,imu,inu)
          end do
        end do

        ! Rotate dipole vector to finally obtain dipole-vector in crystal coordinates 
        !(2-nd rotation; the orbitals are already rotated in "doscentros")
        do ix = 1, 3
          do iy = 1, 3
            dipc(ix,:,:,ineigh,iatom) = dipc(ix,:,:,ineigh,iatom) + eps(ix,iy) *dipcm(iy,:,:)
          enddo
        enddo

        !  DERIVATIVES: Compute dippc from dippcm, eps and deps. First index (ix): x,y,z:
        !  derivation coordinate. Second index: component of the dipole (X,Y,Z). third
        !  e.g. dippc(2,3,:,:,:,:) means the derivative with respect to y of the z-component of
        !  the dipole vector;
        !  fourth indexes: mu,nu orbitals.  
        !  The final value is dippc.
        !  dippcm from "doscentros" is an intermediate value in which the orbitals are rotated 
        !  but not the dipole-vector 
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
