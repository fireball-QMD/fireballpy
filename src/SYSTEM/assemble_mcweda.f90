! ===========================================================================
!       This routine reads the different data file.
! ===========================================================================
subroutine assemble_mcweda ()
  use M_system
  implicit none
  integer iatom
  integer jatom
  integer ineigh
  integer mbeta
  integer kforce

  if (Kscf .eq. 1) then
    call initneighbors (natoms, ivdw, nstepi)
    call num_neigh_tot (numorb_max)
    call backnay ()
    call neighbors_pairs(icluster)
    call common_neighbors (nprocs, my_proc, iordern, iwrtneigh_com)
    call common_neighborsPP (nprocs, my_proc, iordern, iwrtneigh_com, icluster)
  end if ! end if (Kscf .eq. 1)

  ! ewald energy
  kforce = 0
  call get_ewald (nprocs, my_proc, kforce, icluster, itheory, iordern)

  ! A S S E M B L E    H A M I L T O N I A N
  ! O B T A I N   B A N D - S T R U C T U R E
  !-------------------------------------------
  ! Assemble the matrix elements of the Hamiltonian - all in real space.
  ! Set up neigh_self.  The variable neigh_self(natoms) is the ineigh value
  ! for the "self interaction".  Find the neighbor-number of iatom with itself
  ! (neigh_self) in order to put the result of VNA_atom (doscentros) into
  ! VNA(mu,nu,iatom,neigh_self).
  ! Initialize to something ridiculous.

  neigh_self = -999
  do iatom = 1, natoms
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) neigh_self(iatom) = ineigh
    end do
  end do
  neighPP_self = -999
  do iatom = 1, natoms
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0)  neighPP_self(iatom) = ineigh
   end do
  end do

  ! assemble_1c
  ! Assemble the one-center exchange-correlation interactions.
  call assemble_olsxc_1c (natoms, itheory, iforce)

  if (V_intra_dip .eq. 1) call assemble_1c_vdip (iforce)

  !  ----------------- assemble_2c ------------------------
  ! We now do the 2center contributions - assemble_2c does the assembly and
  ! doscentros does the interpolating and "fundamental" calculations.
  ! assemble_2c ONLY does 2c terms. No 3c terms allowed. See assemble_3c
  ! and trescentros for 3c terms.
  if (Kscf .eq. 1) then
    call assemble_sVNL (iforce)
    call assemble_2c (iforce)
    !AQUI
    call assemble_2c_PP (nprocs, iforce, iordern)
  end if ! end if of Kscf = 1

! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).
          if (itheory_xc .eq. 1 ) then
          !write (*,*) ' Assemble SN-xc exchange-correlation interactions. '

          if (itheory .eq. 1) then
           call average_ca_rho (nprocs, Kscf, iforce, iordern, igauss)
          else
           call average_rho (nprocs, Kscf, iforce, iordern, igauss)
          endif

           !write (*,*) ' Assembling on-site part.'
           call assemble_snxc_on (natoms, nprocs, my_proc, iordern, itheory, &
     &                            uxcdcc_sn)

           !write (*,*) ' Assembling off-site part.'
           call assemble_snxc_off (natoms, nprocs, my_proc, iordern,    &
     &                             itheory)
          end if ! if (itheory_xc = 1)

          if (itheory_xc .eq. 2 ) then
           !write (*,*) ' Assemble OLS-xc exchange-correlation interactions.'

           if (itheory .eq. 1) then
            call average_ca_rho (nprocs, Kscf, iforce, iordern, igauss)
            !call average_rho (nprocs, Kscf, iforce, iordern, igauss)
           else
            call average_rho (nprocs, Kscf, iforce, iordern, igauss)
           endif

           !write (*,*) ' Assembling on-site part.'
           call assemble_olsxc_on (natoms, nprocs, my_proc, iordern,    &
     &                             itheory, uxcdcc_ols)

           !write (*,*) ' Assembling off-site part.'
           call assemble_olsxc_off (nprocs, my_proc, iordern, itheory)
          end if ! if (itheory_xc = 2)

!JIMM
          if (itheory .eq. 1) then
           !write (*,*) ' Assemble two-center DOGS interactions. '
           if (idipole .eq. 0) call assemble_ca_2c (nprocs, iforce, iordern)
           if (idipole .eq. 1) call assemble_ca_2c_dip (nprocs, iforce, iordern)
          endif
! ===========================================================================
!                               assemble_3c
! ===========================================================================
! We now do the 3center contributions - assemble_3c does the assembly and
! trecentros does the interpolating and "fundamental" calculations.
! assemble_3c ONLY does 3c terms. No 2c terms allowed. See assemble_2c
! and doscentros for 2c terms.
          !write(*,*) '  '
          if (Kscf .eq. 1) then
           !write (*,*) ' Assemble three-center interactions. '
           call assemble_3c (nprocs, iordern, igauss, itheory_xc)
           !write (*,*) ' Assemble three-center PP interactions. '
           call assemble_3c_PP (nprocs, iordern)
! JIMM
           if (iqmmm .eq.1 ) then
             !write (*,*) ' Assemble qm/mm interactions. '
             if (idipole .eq. 0) call assemble_qmmm (nprocs, iordern)
             if (idipole .eq. 1) call assemble_qmmm_dip (nprocs, iordern)
           else
             eqmmm = 0.0d0
             ewaldqmmm = 0.0d0
           end if
          end if
!JIMM
          if (itheory .eq. 1) then
           !write (*,*) ' Assemble three-center DOGS interactions. '
           if (idipole .eq. 0) call assemble_ca_3c (nprocs, iordern, igauss)
           if (idipole .eq. 1) call assemble_ca_3c_dip (nprocs, iordern, igauss)

! Add assemble_lr here for the long long-range ewald contributions
           !write (*,*) ' Assemble long-range interactions. '
           if (idipole .eq. 0) call assemble_lr (nprocs, iordern)
           if (idipole .eq. 1) call assemble_lr_dip (nprocs, iordern)
          endif

          !write (*,*) ' ***************************************************** '

! ===========================================================================
!                                 Build H
! ===========================================================================
! Set up the full Hamiltonian and !writeout HS.dat.
          call buildh (nprocs, itheory, iordern, itestrange,    &
     &                 testrange, ibias, iwrtHS)
! ===========================================================================
! For iwrthampiece .eq. 1 (file - output.input), !write out Hamiltonian pieces
          if (iwrthampiece .eq. 1) then
           call hampiece (itheory)
          end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine assemble_mcweda

