! This routine calculates the average densities with charge transfer.
subroutine average_ca_rho ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: iforce, xc_overtol, natoms, ratom, imass, neigh_max, Kscf, neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, &
    & neigh_comm, neigh_comn, neigh_back, numorb_max, Qin, rho_off, rhoij_off, sm_mat, spm_mat, rho_on, arho_on, rhoi_on, &
    & arhoi_on, arhop_on, rhop_on, arhoij_off, arho_off, arhopij_off, arhop_off, rhop_off, rhopij_off, xl
  use M_fdata, only: nssh, num_orb, nsh_max
  implicit none
  integer iatom
  integer ibeta
  integer imu
  integer in1
  integer in2
  integer in3
  integer indna
  integer ineigh
  integer interaction
  integer interaction0
  integer inu
  integer isorp, ialp
  integer issh
  integer jatom
  integer jneigh
  integer jbeta
  integer jssh
  integer mbeta
  integer mneigh
  real(double) cost
  real(double) x
  real(double) y
  real(double), dimension (3, 3, 3) :: deps
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r21
  real(double), dimension (3) :: rhat
  ! crystaline matrices
  real(double), dimension (numorb_max, numorb_max) :: rho_2c
  real(double), dimension (numorb_max, numorb_max) :: rhoi_2c
  ! molecular spehrical matrices
  real(double), dimension (numorb_max, numorb_max) :: rhomx
  real(double), dimension (3, numorb_max, numorb_max) :: rhompx
  real(double), dimension (nsh_max, nsh_max) :: rhom_2c
  real(double), dimension (nsh_max, nsh_max) :: rhomi_2c
  real(double), dimension (3, nsh_max, nsh_max) :: rhomp_2c
  real(double), dimension (:, :, :, :), allocatable :: rhom_3c
  real(double), dimension (nsh_max, nsh_max) :: rhomm
  real(double), dimension (3, nsh_max, nsh_max) :: rhompm
  real(double), dimension (3) :: rnabc
  real(double), dimension (3) :: rna
  real(double), dimension (3) :: sighat
  real(double), dimension (nsh_max, nsh_max) :: sm
  real(double), dimension (3, nsh_max, nsh_max) :: spm
  !
  real(double) rho_modified

  rhomx=0.0d0
  rhompx=0.0d0

  allocate (rhom_3c (nsh_max, nsh_max, neigh_max, natoms))

  if (Kscf .eq. 1) sm_mat = 0.0d0
  if (Kscf .eq. 1 .and. iforce .eq. 1) spm_mat = 0.0d0
  
  !   -----  ON SITE PART  ------
  ! We assemble on-site density matrices
  ! <mu i| rho_i^o | nu i>  ......  rhoi0_on; arhoi0_on
  ! <mu i| rho_j | nu i>    ......  rho_on; arho_on

  arho_on   = 0.0d0
  arhoi_on  = 0.0d0
  rho_on    = 0.0d0
  rhoi_on   = 0.0d0
  ! forces
  arhop_on = 0.0d0
  rhop_on  = 0.0d0

  do iatom = 1,natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)

    rho_2c   = 0.0d0
    rhoi_2c  = 0.0d0
    rhom_2c   = 0.0d0
    rhomi_2c  = 0.0d0
    ! spherical overlap <i,mu| i,nu>
    y = 0.0d0
    in2 = in1
    isorp = 0
    interaction0 = 23
    in3 = in2
    sighat = 0.0d0
    eps = 0.0d0

    do issh = 1,3
      eps(issh,issh) = 1.0d0
    enddo
    call doscentrosS (interaction0, isorp, iforce, in1, in2, in3, y, eps, sm, spm)

    do ineigh = 1, neighn(iatom) 
      rhomp_2c  = 0.0d0
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

      ! CALL DOSCENTROS AND GET VXC FOR ATM CASE - AVERAGE DENSITY APPROXIMATION
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        interaction = 17
        interaction0 = 22
        in3 = in1
        do isorp = 1, nssh(in2)
          call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
          call doscentrosS (interaction0, isorp, iforce, in1, in2, in3, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in3)
              ! scf atomic density term
              rhoi_on(imu,inu,iatom) = rhoi_on(imu,inu,iatom) + rhomx(imu,inu)*Qin(isorp,jatom)
              rhoi_2c(imu,inu) = rhoi_2c(imu,inu) + rhomx(imu,inu)*Qin(isorp,jatom)
              ! scf onsite density term
              rho_on(imu,inu,iatom) = rho_on(imu,inu,iatom) + rhomx(imu,inu)*Qin(isorp,jatom)
              rho_2c(imu,inu) = rho_2c(imu,inu) + rhomx(imu,inu)*Qin(isorp,jatom)
            end do
          end do
          do inu = 1, nssh(in1)
            do imu = 1, nssh(in3)
              rhom_2c(imu,inu) = rhom_2c(imu,inu) + rhomm(imu,inu)*Qin(isorp,jatom)
              rhomi_2c(imu,inu) = rhomi_2c(imu,inu) + rhomm(imu,inu)*Qin(isorp,jatom)
            end do   ! endo imu
          end do   ! enddo inu
        end do   ! endo do isorp
      else
        interaction = 17
        interaction0 = 22
        in3 = in1
        do isorp = 1, nssh(in2)
          call doscentros (interaction, isorp, iforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
          call doscentrosS (interaction0, isorp, iforce, in1, in2, in3, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in3)
              rho_on(imu,inu,iatom) = rho_on(imu,inu,iatom) + rhomx(imu,inu)*Qin(isorp,jatom)
              rho_2c(imu,inu) = rho_2c(imu,inu) + rhomx(imu,inu)*Qin(isorp,jatom)
              rhop_on(:,imu,inu,ineigh,iatom) = rhop_on(:,imu,inu,ineigh,iatom) + rhompx(:,imu,inu)*Qin(isorp,jatom)
            end do 
          end do 
          ! spherical symetric
          do inu = 1, nssh(in1)
            do imu = 1, nssh(in3)
              rhom_2c(imu,inu) = rhom_2c(imu,inu) + rhomm(imu,inu)*Qin(isorp,jatom)
              rhomp_2c(:,imu,inu) = rhomp_2c(:,imu,inu) +  rhompm(:,imu,inu)*Qin(isorp,jatom)
            end do   
          end do   
        end do   !isorp
      end if   ! end if (iatom.eq.jatom)
      ! Now assemble the derivative average density using the density pieces from above.
      do issh = 1,nssh(in1)
        do jssh = 1,nssh(in1)
          if(abs(sm(issh,jssh)) .lt. xc_overtol) then
            if (sm(issh,jssh) .gt. 0.0d0) then
              sm(issh,jssh) =  xc_overtol
            else
              sm(issh,jssh) =  -1.0d0*xc_overtol
            endif
          endif
          arhop_on(:,issh,jssh,ineigh,iatom) =  arhop_on(:,issh,jssh,ineigh,iatom) &
          &     + ( rhomp_2c(:,issh,jssh)*sm(issh,jssh) - rhom_2c(issh,jssh)*spm(:,issh,jssh) ) / ( sm(issh,jssh)*sm(issh,jssh) )
        enddo   ! do jssh
      enddo   ! do issh
    end do   ! end do ineigh

    ! We assemble the average density in molecular coordinates.  The average density is used later in the exchange-correlation energy, potential, corresponding derivatives.
    ! Now assemble the average density using the density pieces from above. Loop over shells.
    do issh = 1,nssh(in1)
      do jssh = 1,nssh(in1)
        if(abs(sm(issh,jssh)) .lt. xc_overtol) then
          if (sm(issh,jssh) .gt. 0.0d0) then
            sm(issh,jssh) =  xc_overtol
          else
            sm(issh,jssh) =  -1.0d0*xc_overtol
          endif
        endif
        arho_on(issh,jssh,iatom) = (arho_on(issh,jssh,iatom)  + rhom_2c(issh,jssh) / sm(issh,jssh))
        arhoi_on(issh,jssh,iatom) = (arhoi_on(issh,jssh,iatom) + rhomi_2c(issh,jssh) / sm(issh,jssh))
      enddo   ! do jssh
    enddo   ! do issh
  enddo   ! end do iatom
  !   -----  OFF SITE PART  ------
  ! We assemble off-site density matrices
  ! <mu i| (rho_i+rho_j) | nu j>  ......  rhoij_off; arhoij_off
  ! <mu i| rho_k | nu j>   ......  rho_off; arho_off
  ! Initialize arrays to zero.
  arhoij_off  = 0.0d0
  arhopij_off = 0.0d0
  arho_off   = 0.0d0
  arhop_off  = 0.0d0
  rho_off = 0.0d0
  rhop_off = 0.0d0
  rhoij_off = 0.0d0
  rhopij_off = 0.0d0
  rhom_3c = 0.0d0

  !    T H R E E - C E N T E R   P A R T
  do ialp = 1, natoms
    rna(:) = ratom(:,ialp)
    indna = imass(ialp)
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh,ialp)
      if (mneigh .ne. 0) then
        iatom = neigh_comj(1,ineigh,ialp)
        ibeta = neigh_comb(1,ineigh,ialp)
        r1(:) = ratom(:,iatom) + xl(:,ibeta)
        in1 = imass(iatom)
        jatom = neigh_comj(2,ineigh,ialp)
        jbeta = neigh_comb(2,ineigh,ialp)
        r2(:) = ratom(:,jatom) + xl(:,jbeta)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,mneigh)
        r21(:) = r2(:) - r1(:)
        y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
        if (y .lt. 1.0d-05) then
          sighat(1) = 0.0d0
          sighat(2) = 0.0d0
          sighat(3) = 1.0d0
        else
          sighat(:) = r21(:)/y
        end if
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        if (x .lt. 1.0d-05) then
          rhat(1) = 0.0d0
          rhat(2) = 0.0d0
          rhat(3) = 0.0d0
        else
          rhat(:) = rnabc(:)/x
        end if
        cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
        call epsilon (rhat, sighat, eps)
        interaction = 3
        do isorp = 1, nssh(indna)
          call trescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, rhomx)
          call trescentrosS (isorp, in1, in2, indna, x, y, cost, rhomm)
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              rho_off(imu,inu,mneigh,iatom) = rho_off(imu,inu,mneigh,iatom) + rhomx(imu,inu)*Qin(isorp,ialp)
              rho_off(inu,imu,jneigh,jatom) = rho_off(imu,inu,mneigh,iatom)
            end do
          end do
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              rhom_3c(imu,inu,mneigh,iatom) = rhom_3c(imu,inu,mneigh,iatom) + rhomm(imu,inu)*Qin(isorp,ialp)
              rhom_3c(inu,imu,jneigh,jatom) = rhom_3c(imu,inu,mneigh,iatom) 
            end do
          end do
        end do   ! do isorp
      end if
    end do   ! end do neigh_comn(ialp)
  end do    ! end do ialp
  !      T W O - C E N T E R   P A R T
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
         ! Do nothing here - special case. Interaction already calculated in on-site
         ! stuff, i.e. assemble_olsxc_on.f90 We calculate only off diagonal elements here.
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
        interaction = 15
        interaction0 = 20
        in3 = in1
        rhom_2c = 0.0d0
        rhomp_2c = 0.0d0
        do isorp = 1, nssh(in3)
          call doscentros (interaction, isorp, iforce, in1, in3, in2, y, eps, deps, rhomx, rhompx)
          call doscentrosS (interaction0, isorp, iforce, in1, in3, in2, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              rhoij_off(imu,inu,ineigh,iatom) = rhoij_off(imu,inu,ineigh,iatom)        + rhomx(imu,inu)*Qin(isorp,iatom)
              rhopij_off(:,imu,inu,ineigh,iatom) =  rhopij_off(:,imu,inu,ineigh,iatom) + rhompx(:,imu,inu)*Qin(isorp,iatom)
              rho_off(imu,inu,ineigh,iatom)      = rho_off(imu,inu,ineigh,iatom)       + rhomx(imu,inu)*Qin(isorp,iatom)
              rhop_off(:,imu,inu,ineigh,iatom)   =  rhop_off(:,imu,inu,ineigh,iatom)   + rhompx(:,imu,inu)*Qin(isorp,iatom)
            end do
          end do
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              rhom_2c(imu,inu) = rhom_2c(imu,inu) + rhomm(imu,inu)*Qin(isorp,iatom)
              rhomp_2c(:,imu,inu) = rhomp_2c(:,imu,inu) + rhompm(:,imu,inu)*Qin(isorp,iatom)
            end do
          end do
        end do

        interaction = 16
        interaction0 = 21
        in3 = in2
        do isorp = 1, nssh(in3)
          call doscentros (interaction, isorp, iforce, in1, in3, in2, y, eps, deps, rhomx, rhompx)
          call doscentrosS (interaction0, isorp, iforce, in1, in3, in2, y, eps, rhomm, rhompm)
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              rhoij_off(imu,inu,ineigh,iatom) = rhoij_off(imu,inu,ineigh,iatom) + rhomx(imu,inu)*Qin(isorp,jatom)
              rhopij_off(:,imu,inu,ineigh,iatom) =  rhopij_off(:,imu,inu,ineigh,iatom)  + rhompx(:,imu,inu)*Qin(isorp,jatom)
              rho_off(imu,inu,ineigh,iatom) =  rho_off(imu,inu,ineigh,iatom) + rhomx(imu,inu)*Qin(isorp,jatom)
              rhop_off(:,imu,inu,ineigh,iatom) = rhop_off(:,imu,inu,ineigh,iatom) + rhompx(:,imu,inu)*Qin(isorp,jatom)
            end do
          end do
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              rhom_2c(imu,inu) = rhom_2c(imu,inu) + rhomm(imu,inu)*Qin(isorp,jatom)
              rhomp_2c(:,imu,inu) =  rhomp_2c(:,imu,inu) + rhompm(:,imu,inu)*Qin(isorp,jatom)
            end do
          end do
        end do
        if (Kscf .eq. 1) then
          isorp = 0
          interaction0 = 23
          in3 = in2
          call doscentrosS (interaction0, isorp, iforce, in1, in2, in3, y, eps, sm, spm)
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              sm_mat(imu,inu,ineigh,iatom) = sm(imu,inu)
              if (iforce .eq. 1) spm_mat(:,imu,inu,ineigh,iatom) = spm(:,imu,inu)
            end do
          end do
        else
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              sm(imu, inu) = sm_mat(imu, inu, ineigh, iatom)
              spm(:, imu, inu) = spm_mat(:, imu, inu, ineigh, iatom)
            end do
          end do
        end if
        rho_modified= 0.d0
     
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            if (abs(sm(issh,jssh)) .lt. xc_overtol) then
              if (sm(issh,jssh) .gt. 0.0d0) then
                sm(issh,jssh) =  xc_overtol
              else
                sm(issh,jssh) =  -1.0d0*xc_overtol
              endif
            endif
            arhoij_off(issh,jssh,ineigh,iatom) =  arhoij_off(issh,jssh,ineigh,iatom) +  rhom_2c(issh,jssh)/sm(issh,jssh)
            arhopij_off(:,issh,jssh,ineigh,iatom) = arhopij_off(:,issh,jssh,ineigh,iatom)  + (rhomp_2c(:,issh,jssh)*sm(issh,jssh) - rhom_2c(issh,jssh)*spm(:,issh,jssh)) /(sm(issh,jssh)*sm(issh,jssh))
            rho_modified=rhom_3c(issh,jssh,ineigh,iatom)/sm(issh,jssh)
            arho_off(issh,jssh,ineigh,iatom) =  arho_off(issh,jssh,ineigh,iatom)  + rhom_2c(issh,jssh)/sm(issh,jssh) + rho_modified
            arhop_off(:,issh,jssh,ineigh,iatom) = arhop_off(:,issh,jssh,ineigh,iatom) + (rhomp_2c(:,issh,jssh)*sm(issh,jssh) - rhom_2c(issh,jssh)*spm(:,issh,jssh)) /(sm(issh,jssh)*sm(issh,jssh))
          end do   
        end do   
      end if
    end do   ! do ineigh
  end do    ! do iatom

  deallocate (rhom_3c)
end subroutine average_ca_rho
