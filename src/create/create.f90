! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.

! create.f
! Program Description
! ======================================================================
!       Create matrix element interactions as a function of atomtype and
! geometry
!
! ======================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
! and
! Cyra Roldán Piñero
! Departamento de Física Teórica de la Materia Condensada and
!   Condensed Matter Institute
! Universidad Autónoma de Madrid

! ======================================================================
!
! Program Declaration
! ======================================================================
program create
  use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit, stdout => output_unit
  use :: xc, only: xc_init, xc_end
  use :: wavefunctions, only: wf_init, wf_end
  use :: potentials, only: pot_init, pot_end
  use :: pseudopotentials, only: pp_init, pp_end
  use :: onecenter, only: onecenter_calc, ONECENTER_NPOINTS, &
    &                     ONECENTER_XC, &
    &                     ONECENTER_GOVERLAP1, &
    &                     ONECENTER_GOVERLAP2
  use :: twocenter, only: twocenter_calc, TWOCENTER_NPOINTS_D, TWOCENTER_NPOINTS_Z, TWOCENTER_NPOINTS_RHO, &
    &                     TWOCENTER_DENS_ATOM, &
    &                     TWOCENTER_DENS_ONTOP, &
    &                     TWOCENTER_OVERLAP, &
    &                     TWOCENTER_VNA_ATOM, &
    &                     TWOCENTER_VNA_ONTOP, &
    &                     TWOCENTER_VPP, &
    &                     TWOCENTER_VXC, &
    &                     TWOCENTER_DIP_Z, &
    &                     TWOCENTER_DIP_X, &
    &                     TWOCENTER_DIP_Y, &
    &                     TWOCENTER_COULOMB, &
    &                     TWOCENTER_DENS_ATOM_SPH, &
    &                     TWOCENTER_DENS_ONTOP_SPH, &
    &                     TWOCENTER_OVERLAP_SPH
  implicit none

  include 'parameters.inc'
  include 'exchange.inc'

! The file "quadrature.inc" specifies all the quadrature and grid
! parameters for all of the matrix elements.
  include 'quadrature.inc'
  include 'wavefunctions.inc'
  include 'pseudopotentials.inc'
! Argument Declaration and Description
! ======================================================================

! Local Parameters and Data Declaration
! ======================================================================
  real(kind=dp) abohr
  parameter (abohr = 0.529177249_dp)

! Local Variable Declaration and Description
! ======================================================================
  integer ideriv
  integer iderivmin
  integer iderivmax
  integer iderivtmp
  integer iexc
  integer iexc_new
  integer index
  integer index_max
  integer index_maxPP
  integer index_maxsp
  integer interaction
  integer ispec
  integer isorp
  integer isorpmax
  integer isorpmin
  integer ispnum
  integer istyle
  integer issh
  integer itmp
  integer itype1, itype2, itype3
  integer nstyles

! Needed for determing expansion of Legendre polynomials.
  real(kind=dp) ctheta (ntheta_max)
  real(kind=dp) ctheta_weights (ntheta_max)

! The two center interactions are defined as follows:
! ontop => orbitals at two different sites
! atom-atom => orbitals at same site
! interaction = 0 => do average density
! interaction = 1 => do overlap
! interaction = 2 => do neutral atom/ontop function => vnnaofr
! interaction = 3 => do neutral atom/atom
! interaction = 4 => do non-local
! interaction = 5 => do xc ontop function => vxc
! interaction = 6 => do xc atom function => dvxc
! interaction = 7 => do xc double counting correction function => dexc
! interaction = 8 => do z-dipole
! interaction = 9 => do y-dipole
! interaction = 10 => do x-dipole
! interaction = 11 => do coulomb
! switch.input variables
  integer ibcna
  integer ibcxc
  integer ikinetic
  integer iswitch (0:11)     ! two-center interactions
  integer imuxc1c
  integer inuxc1c
  integer inuxc2c
  integer isnuxc1c
  integer isnuxc2c
  integer V_intra_dip

! theory.input variables
  integer idogs
  integer igauss
  integer iharris
  integer ihubbard
  integer ioomethod
  integer ispin
  integer itest
  integer ngauss
  integer ixc_opt   ! jel-oslxc
  logical ispher

  integer looper2
  integer looper23
  integer looper24
  integer looper25
  integer looper3
  integer looper3a
  integer minderiv2c
  integer maxderiv2c
  integer nspec
  integer nssh2
  integer nssh2PP
  integer nzx1, nzx2
  integer itype

!
! information relevant to the storage of the non-zero matrix elements
! for two and three center matrix elements.
!
! terms 1 .. index_max2c(in1,in2): all two-center matrix elements
!
! terms index_max2c(in1,in2) + 1 .. index_max3c(in1,in2):
!       all additional terms appearing in the three center case
!
  integer index_max2c (nspec_max, nspec_max)
  integer index_max2cPP (nspec_max, nspec_max)
  integer index_max3c (nspec_max, nspec_max)
  integer lleft (nspec_max, nspec_max, inter_max)
  integer lleftPP (nspec_max, nspec_max, inter_max)
  integer lright (nspec_max, nspec_max, inter_max)
  integer lrightPP (nspec_max, nspec_max, inter_max)
  integer mleft (nspec_max, nspec_max, inter_max)
  integer mleftPP (nspec_max, nspec_max, inter_max)
  integer mright (nspec_max, nspec_max, inter_max)
  integer mrightPP (nspec_max, nspec_max, inter_max)
  integer nleft (nspec_max, nspec_max, inter_max)
  integer nleftPP (nspec_max, nspec_max, inter_max)
  integer nright (nspec_max, nspec_max, inter_max)
  integer nrightPP (nspec_max, nspec_max, inter_max)

! general case
  integer l1 (inter_max)
  integer l1PP (inter_max)
  integer l2 (inter_max)
  integer l2PP (inter_max)
  integer m1 (inter_max)
  integer m1PP (inter_max)
  integer m2 (inter_max)
  integer m2PP (inter_max)
  integer n1 (inter_max)
  integer n1PP (inter_max)
  integer n2 (inter_max)
  integer n2PP (inter_max)
  integer n1sp (1), l1sp (1), m1sp (1)
  integer n2sp (1), l2sp (1), m2sp (1)
  integer lalpha
  integer malpha
  integer nalpha

  integer lssh (nspec_max, nsh_max)      ! l quantum number
  integer lsshPP (nspec_max, nsh_max)    ! l quantum number for PP
  integer nssh (nspec_max)               ! number of shells
  integer nsshPP (nspec_max)             ! number of shells for PP
  integer nzx (nspec_max)

  real(kind=dp) fraction
  real(kind=dp) max_diff
  real(kind=dp) rcutoff1
  real(kind=dp) rcutoff2
  real(kind=dp) rcutoff3

! distances needed for three center integrals
  real(kind=dp) dbc                      ! maximal bond charge distance
  real(kind=dp) dna                      ! maximal neutral atom distance

  real(kind=dp) etotatom (nspec_max)
  real(kind=dp) rcutoff (nspec_max, nsh_max) ! cutoff radius in bohr
  real(kind=dp) rcutoffa (nspec_max, nsh_max)! cutoff radius in angstroms
  real(kind=dp) rcutoffa_max (nspec_max)     ! cutoff radius in angstroms
  real(kind=dp) xmass (nspec_max)

  character(len=2) atom1, atom2, atom3
  character(len=70) signature
  character(len=45) time_message
  character(len=70) what1, what2, what3

  character(len=2) atom (nspec_max)
  character(len=25) napot (nspec_max, 0:nsh_max)
  character(len=25) ppfile (nspec_max)
  character(len=25) wavefxn (nspec_max, nsh_max)
  character(len=70) what (nspec_max)

  logical read_info

! MPI
  logical iammaster, iammpi
  integer my_proc, nproc

! Procedure
! ======================================================================
! Setup the clm coefficients
  call setup_clm

! Initialize MPI
  call init_MPI (iammaster, iammpi, my_proc, nproc)

  if (iammaster) then
    write (*,*) '  '
    write (*,100)
    write (*,*) '  '
    write (*,*) ' ****** Welcome to the creator package ****** '
    write (*,*) '          running on ',nproc,' nodes'
    write (*,100)
    write (*,*) '  '
    write (*,*) '               Fireballs-2005 '
    write (*,*) '      A fast local orbital QMD Package '
    write (*,*) '  '
    write (*,100)
    write (*,*) '  '

! Identification of the user:
!         write (*,*) '  '
!         write (*,*) ' Please insert your name and other messages. '
!         write (*,*) '  '
    signature = ' fireballpy '
  end if ! end master

!  MPI Broadcast signature
  call signature_MPI (signature)

! **********************************************************************
! Read in create.input and set up some data
  call readcreate (nspec, iammaster, iammpi, atom, what, nssh, &
  &                   lssh, nzx, rcutoff, rcutoffa, rcutoffa_max, &
  &                   xmass, ppfile, napot, wavefxn)

! **********************************************************************
! We now read in a theory.input file. This determines certain defaults
! for the switches which is dependent upon the the level of theory that
! is chosen.
  call readtheory (iammaster, ibcna, ibcxc, ikinetic, iswitch,  &
  &                   imuxc1c, inuxc1c, inuxc2c, isnuxc1c, &
  &                   isnuxc2c, V_intra_dip,  &
  &                   itest, idogs, iharris, ihubbard, ispin,  &
  &                   ioomethod, ixc_opt, igauss, ngauss)

! **********************************************************************

  if (iammaster) write (*,*) ' Reading wavefunctions. '
  do ispec = 1, nspec
    do issh = 1, nssh(ispec)
      call readpsi (ispec, issh, lssh(ispec,issh), &
      &                  rcutoff(ispec,issh), xnocc(issh,ispec), &
      &                  nzx(ispec), wavefxn(ispec,issh), iammaster)
    end do
  end do

! Read in the pseudopotential
  if (iammaster) write (*,*) ' Reading pseudopotentials. '
  do ispec = 1, nspec
! jel-PP
    call readvpp (ispec, ppfile(ispec), nsshPP, lsshPP, iexc, &
    &                 fraction, iammaster)
  end do

  if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
  &      .or. iexc .eq. 10) then
    stop 1
  end if
  if (iexc .ne. 11 .and. (ispin .eq. 1 .or. isnuxc1c .eq. 1 .or. &
  &                                            isnuxc2c .eq. 1)) then
    if (iammaster) then
      write (*,*) ' The exchange-correlation option that you chose '
      write (*,*) ' has not been implemented for the '
      write (*,*) ' spin-polarization case. Only iexc = 11 works! '
      write (*,*) ' Choose iexc = 11 or ispin = 0, and restart. '
    end if
    stop 1
  end if

! Read the neutral and non-neutral atom potentials
  if (iammaster) write (*,*) ' Reading potentials.'
  do ispec = 1, nspec
    call readvnn (ispec, 0, rcutoffa_max(ispec)/abohr, nzx(ispec), &
    &                 napot(ispec,0), etotatom(ispec), iammaster)
    do issh = 1, nssh(ispec)
      call readvnn (ispec, issh, rcutoff(ispec,issh), &
      &                  nzx(ispec), napot(ispec,issh), etotatom(ispec), &
      &                  iammaster)
    end do
  end do

! **********************************************************************
! Call the following subroutines to determine the needed matrix elements
! for each itype1, itype2 pair.
  do itype1 = 1, nspec
    do itype2 = 1, nspec
      call mk_index (itype1, itype2, nspec_max, nsh_max, inter_max, &
      &                   nssh, lssh, nleft, lleft, mleft, nright, &
      &                   lright, mright, index_max2c, index_max3c)

! Write out the results of mk_index
      if (index_max3c(itype1,itype2) .gt. inter_max) then
        write (*,*) ' index_max3c(itype1,itype2) = ',  &
        &                   index_max3c(itype1,itype2)
        write (*,*) ' inter_max = ', inter_max
        write (*,*) ' Redimension index_max in parameters.inc! '
        stop 1
      end if

! For the pseudopotential call a different mk_index routine
      call mk_indexPP (itype1, itype2, nspec_max, nsh_max, &
      &                     inter_max, nssh, nsshPP, lssh, lsshPP, &
      &                     nleftPP, lleftPP, mleftPP, nrightPP, &
      &                     lrightPP, mrightPP, index_max2cPP)
    end do
  end do

! ****************************************************************************
! Write out the info file.
! Only do this on the master
  if (iammaster) then
    inquire (file = 'coutput/info.dat', exist = read_info)
    open (unit = 12, file = 'coutput/info.dat', status = 'unknown')

    time_message = ' No time available! '
    open (unit = 45, file = 'timecreate', status = 'unknown')
    do itmp = 1, 10000
      read (45, 106, err = 445, end = 445) time_message
    end do
445 continue
    close (unit = 45)

    if (.not. read_info) then
      write (12,104) signature
      write (12,*) nspec, ' - Number of species '

      do ispec = 1, nspec
        write (12,100)
        write (12,301) ispec
        write (12,302) atom(ispec)
        write (12,303) nzx(ispec)
        write (12,304) xmass(ispec)
        write (12,305) nssh(ispec)
        write (12,306) (lssh(ispec,issh), issh = 1, nssh(ispec))
        write (12,307) nsshPP(ispec)
        write (12,308) (lsshPP(ispec,issh), issh = 1, nsshPP(ispec))
!jel-PP
        write (12,314) rcPP(ispec)
        write (12,309) (xnocc(issh,ispec), issh = 1, nssh(ispec))
        write (12,310) (rcutoff(ispec,issh), issh = 1, nssh(ispec))
        write (12,311) (wavefxn(ispec,issh), issh = 1, nssh(ispec))
        write (12,312) (napot(ispec,issh), issh = 0, nssh(ispec))
        write (12,313) etotatom(ispec)
        write (12,100)
      end do

! If the file already exists then just write out the time
    else
      do itmp = 1, 10000
        read (12, *, err = 446, end = 446)
      end do
446   continue
      backspace 12
    end if

    write (12, 106) time_message
    write (12, 107) itest, iharris, idogs, ihubbard, ispin, &
    &                   ioomethod, ixc_opt
    write (12, 108) igauss, ngauss
    write (12, 109) imuxc1c, ikinetic, &
    &                   (iswitch(interaction), interaction = 0, 11), &
    &                   ibcna, ibcxc, inuxc1c, inuxc2c, isnuxc1c, &
    &                   isnuxc2c
    close (unit = 12)

! Done with setup, now get to work
! ======================================================================

    ! TODO: be able to pick whatever iexc we want
    if (iexc == 3) then
      call xc_init(iexc1=1, iexc2=9)
    else if (iexc == 9) then
      call xc_init(iexc1=106, iexc2=131)
    else
      write (stderr, "(a)") "TODO: be able to pick whatever iexc we want"
      error stop 1
    end if
    call wf_init()
    call pot_init()
    call pp_init()


    write (stdout, "(a)") "  ==== ONE CENTER INTEGRALS ===="
    write (stdout, "(2x,a,i6)") "Number of rho integration points:", ONECENTER_NPOINTS
    write (stdout, "(a)") repeat("=", 30)
    ! TODO: mpi
    do itype1 = 1, nspec
      call onecenter_calc(ONECENTER_XC + &
        &                 ONECENTER_GOVERLAP1 + &
        &                 ONECENTER_GOVERLAP2, ispec=itype1)
    end do
    write (stdout, "(2x,a)") repeat("=", 30)
    write (stdout, "(a)") ""

! JPL 1999 Exact exchange interactions.
    ! if (imuxc1c .eq. 1 .and. iexc .eq. 12) then
    !   call x_1c (nsh_max, nspec, nspec_max, fraction, nssh, lssh, &
    !   &               drr_rho, rcutoffa_max, what, signature)
    ! end if

    write (stdout, "(a)") "  ==== TWO CENTER INTEGRALS ===="
    write (stdout, "(2x,a,i6)") "Number of grid points:", TWOCENTER_NPOINTS_D
    write (stdout, "(2x,a,i6)") "Number of z integration points:", TWOCENTER_NPOINTS_Z
    write (stdout, "(2x,a,i6)") "Number of rho integration points:", TWOCENTER_NPOINTS_RHO
    write (stdout, "(2x,a)") repeat("=", 30)
    ! TODO: mpi
    do itype1 = 1, nspec
      do itype2 = 1, nspec
        call twocenter_calc(TWOCENTER_DENS_ATOM + &
          &                 TWOCENTER_DENS_ONTOP + &
          &                 TWOCENTER_OVERLAP + &
          &                 TWOCENTER_VNA_ATOM + &
          &                 TWOCENTER_VNA_ONTOP + &
          &                 TWOCENTER_VPP + &
          &                 TWOCENTER_VXC + &
          &                 TWOCENTER_DIP_Z + &
          &                 TWOCENTER_DIP_X + &
          &                 TWOCENTER_DIP_Y + &
          &                 TWOCENTER_COULOMB + &
          &                 TWOCENTER_DENS_ATOM_SPH + &
          &                 TWOCENTER_DENS_ONTOP_SPH + &
          &                 TWOCENTER_OVERLAP_SPH, ispec=itype1, jspec=itype2)
      end do
    end do
    write (stdout, "(a)") repeat("=", 30)
    write (stdout, "(a)") ""

  end if ! end master

! ======================================================================
! I. Perform three-center calculations
! ======================================================================
! Note igauss controls whether or not we ACTUALLY compute the gaussian
! integrals EVEN WHEN REQUESTED BY igauss3C. The point is that switch
! allows you to NOT compute things EVEN IF YOU NEED THEM!
!        if (igauss .eq. 1) call gausscreate (ngauss)

! We have calculated BOTH 3XC and 3NA gaussian fits. So we need not
! bother with the 3XC or the 3NA parts below.
! We are not doing gaussian fits of anything.

! **********************************************************************
!
!  =====>         1b. Three center exchange-correlation matrix element
!                                   (three-center SNXC)
!
! **********************************************************************

  ideriv = 0

  if (ibcxc .eq. 1 .and. ixc_opt .eq. 1 ) then
    if (iammaster) then
      write (*,*) ' Calculating three-center exchange-correlation '
      write (*,*) ' interactions (SNXC). '
    end if
    nstyles = ideriv
    call gleg (ctheta, ctheta_weights, ntheta_max)
    do looper3a = 1, nspec*nspec*nspec
      if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) &
      &        .eq. my_proc) then
        itmp   = looper3a
        itype3 = 1 + int((itmp - 1)/(nspec*nspec))
        itmp   = itmp - (itype3 - 1)*(nspec*nspec)
        itype2 = 1 + int((itmp - 1)/nspec)
        itmp   = itmp - (itype2 - 1)*nspec
        itype1 = itmp

        rcutoff1 = rcutoffa_max(itype1)
        rcutoff2 = rcutoffa_max(itype2)
        rcutoff3 = rcutoffa_max(itype3)

        atom1 = atom(itype1)
        atom2 = atom(itype2)
        atom3 = atom(itype3)

        what1 = what(itype1)
        what2 = what(itype2)
        what3 = what(itype3)

        dbc = rcutoff1 + rcutoff2
        dna = rcutoff3 + max(rcutoff1,rcutoff2)

        index_max = index_max3c(itype1,itype2)
        do index = 1, index_max
          n1(index) = nleft(itype1,itype2,index)
          l1(index) = lleft(itype1,itype2,index)
          m1(index) = mleft(itype1,itype2,index)
          n2(index) = nright(itype1,itype2,index)
          l2(index) = lright(itype1,itype2,index)
          m2(index) = mright(itype1,itype2,index)
        end do

! If the cutoffs for the atoms are too drastically different, then the grid
! size should be increased. Otherwise, the number of non-zero points may
! be too few.
        max_diff = max(abs(rcutoff1 - rcutoff2), &
        &                    abs(rcutoff2 - rcutoff3), &
        &                    abs(rcutoff3 - rcutoff1))
        if (max_diff .gt. 2.0_dp) then
          write (*,*) ' ************ WARNING ************* '
          write (*,*) ' You have at least two species which have '
          write (*,*) ' rcutoff''s which differ by more than 2.0 '
          write (*,*) ' Angstroms. It is highly advisable that you '
          write (*,*) ' increase the number of mesh points, so as to '
          write (*,*) ' avoid a case where you may end up with many '
          write (*,*) ' zeros, and too few non-zero elements in your '
          write (*,*) ' grid. '
        end if
! xc3c_SN now we call threecenter() with new option 'interaction=3'
        interaction = 3
        ispher = .false.

! Even in case of harris option we need all interactions
        isorpmin = 1
        isorpmax = nssh(itype3)

! Do not parallelize over ispnum, because it involves little work
        ispnum = isorpmax - isorpmin + 1
        call threecenter (itype1, itype2, itype3, index_max, iexc, &
        &                       interaction, nzx, nssh, n1, l1, m1, &
        &                       n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
        &                       atom1, atom2, atom3, what1, what2, what3, &
        &                       dbc, dna, signature, ctheta,  &
        &                       ctheta_weights, isorpmin, &
        &                       isorpmax, ispnum, iammaster, ispher)

      end if ! end MPI which node
    end do
  end if
! xc3c_SN: end of the added part


! **********************************************************************
!
!  =====>         2. Three center neutral atom matrix element
!
! **********************************************************************
  if (ibcna .eq. 1) then
    if (iammaster) then
      write (*,*) ' Calculating three-center neutral atom and '
      write (*,*) ' charged atom interactions. '
    end if
    call gleg (ctheta, ctheta_weights, ntheta_max)
    do looper3a = 1, nspec*nspec*nspec
      if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) &
      &        .eq. my_proc) then
        itmp   = looper3a
        itype3 = 1 + int((itmp - 1)/(nspec*nspec))
        itmp   = itmp - (itype3 - 1)*(nspec*nspec)
        itype2 = 1 + int((itmp - 1)/nspec)
        itmp   = itmp - (itype2 - 1)*nspec
        itype1 = itmp

        rcutoff1 = rcutoffa_max(itype1)
        rcutoff2 = rcutoffa_max(itype2)
        rcutoff3 = rcutoffa_max(itype3)

        atom1 = atom(itype1)
        atom2 = atom(itype2)
        atom3 = atom(itype3)

        what1 = what(itype1)
        what2 = what(itype2)
        what3 = what(itype3)

        dbc = rcutoff1 + rcutoff2
        dna = rcutoff3 + max(rcutoff1,rcutoff2)

        index_max = index_max3c(itype1,itype2)
        do index = 1, index_max
          n1(index) = nleft(itype1,itype2,index)
          l1(index) = lleft(itype1,itype2,index)
          m1(index) = mleft(itype1,itype2,index)
          n2(index) = nright(itype1,itype2,index)
          l2(index) = lright(itype1,itype2,index)
          m2(index) = mright(itype1,itype2,index)
        end do

! If the cutoffs for the atoms are too drastically different, then the grid
! size should be increased. Otherwise, the number of non-zero points may
! be too few.
        max_diff = max(abs(rcutoff1 - rcutoff2), &
        &                    abs(rcutoff2 - rcutoff3), &
        &                    abs(rcutoff3 - rcutoff1))
        if (max_diff .gt. 2.0_dp) then
          write (*,*) ' ************ WARNING ************* '
          write (*,*) ' You have at least two species which have '
          write (*,*) ' rcutoff''s which differ by more than 2.0 '
          write (*,*) ' Angstroms. It is highly advisable that you '
          write (*,*) ' increase the number of mesh points, so as to '
          write (*,*) ' avoid a case where you may end up with many '
          write (*,*) ' zeros, and too few non-zero elements in your '
          write (*,*) ' grid. '
        end if

        interaction = 1
        ispher = .false.
        if (idogs .eq. 1) then
          isorpmin = 1
          isorpmax = nssh(itype3)
        else
          isorpmin = 0
          isorpmax = 0
        end if
        if (iharris .eq. 1) isorpmin = 0

! Do not parallelize over ispnum, because it involves little work
        ispnum = isorpmax - isorpmin + 1
        call threecenter (itype1, itype2, itype3, index_max, iexc, &
        &                       interaction, nzx, nssh, n1, l1, m1, &
        &                       n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
        &                       atom1, atom2, atom3, what1, what2, what3, &
        &                       dbc, dna, signature, &
        &                       ctheta, ctheta_weights, isorpmin, &
        &                       isorpmax, ispnum, iammaster, ispher)
      end if ! end MPI which node
    end do
  end if

! ======================================================================
! II. Compute purely two center cases (overlap, dipole, coulomb
! integrals, and kinetic energy).
! ======================================================================
  do looper2 = 1, nspec*nspec
    if (mod(looper2,nproc) .eq. my_proc) then
      itmp   = looper2
      itype2 = 1 + int((itmp - 1)/nspec)
      itmp   = itmp - (itype2 - 1)*nspec
      itype1 = itmp

      nzx1 = nzx(itype1)
      nzx2 = nzx(itype2)

      rcutoff1 = rcutoffa_max(itype1)
      rcutoff2 = rcutoffa_max(itype2)

      nssh2 = nssh(itype2)

      atom1 = atom(itype1)
      atom2 = atom(itype2)

      what1 = what(itype1)
      what2 = what(itype2)

      index_max = index_max2c(itype1,itype2)
      do index = 1, index_max
        n1(index) = nleft(itype1,itype2,index)
        l1(index) = lleft(itype1,itype2,index)
        m1(index) = mleft(itype1,itype2,index)
        n2(index) = nright(itype1,itype2,index)
        l2(index) = lright(itype1,itype2,index)
        m2(index) = mright(itype1,itype2,index)
      end do

! **********************************************************************
!
!  =====>         1. Kinetic
!
! **********************************************************************
      if (ikinetic .eq. 1) then
        call kinetic (itype1, itype2, atom1, atom2, what1, what2, &
        &                   nzx1, nzx2, rcutoff1, rcutoff2, nssh, lssh, &
        &                   index_max, n1, l1, m1, n2, l2, m2, signature, &
        &                   iammaster)
      end if
    end if
  end do


!
! SPHERIC APPROXIMATION (only OLSXC Method)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  set lssh = lsshxc = 0 ... only s-like orbitals
  lssh = 0
  lsshxc = 0

! **********************************************************************
! Call the following subroutines to determine the needed matrix elements
! for each itype1, itype2 pair.
  do itype1 = 1, nspec
    do itype2 = 1, nspec
      call mk_index (itype1, itype2, nspec_max, nsh_max, inter_max, &
      &                   nssh, lssh, nleft, lleft, mleft, nright, &
      &                   lright, mright, index_max2c, index_max3c)

! Write out the results of mk_index
      if (iammaster) then
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' species 1 = ', itype1, ' species 2 = ', itype2
        write (*,100)
        write (*,*) ' For two-center interactions: '
        write (*,*) ' index_max2c = ', index_max2c(itype1,itype2)
        write (*,500)
        do index = 1, index_max2c(itype1,itype2)
          write (*,501) index, &
          &       nleft(itype1,itype2,index), lleft(itype1,itype2,index), &
          &       mleft(itype1,itype2,index), nright(itype1,itype2,index), &
          &       lright(itype1,itype2,index), mright(itype1,itype2,index)
        end do
        write (*,*) ' Additionally, for three-center interactions: '
        write (*,*) ' index_max3c = ', index_max3c(itype1,itype2)
        write (*,500)
        do index = index_max2c(itype1,itype2) + 1, &
        &                index_max3c(itype1,itype2)
          write (*,501) index, &
          &       nleft(itype1,itype2,index), lleft(itype1,itype2,index), &
          &       mleft(itype1,itype2,index), nright(itype1,itype2,index), &
          &       lright(itype1,itype2,index), mright(itype1,itype2,index)
        end do
      end if ! end master
      if (index_max3c(itype1,itype2) .gt. inter_max) then
        write (*,*) ' index_max3c(itype1,itype2) = ',  &
        &                   index_max3c(itype1,itype2)
        write (*,*) ' inter_max = ', inter_max
        write (*,*) ' Redimension index_max in parameters.inc! '
        stop 1
      end if

    end do
  end do

! **********************************************************************
!
!  =====>         1b. Three center exchange-correlation matrix element
!                                   (three-center OLSXC)
!                                  SPHERICAL APPROXIMATION
! **********************************************************************

  ideriv = 0

  if (ibcxc .eq. 1 .and. ixc_opt .eq. 1 ) then
    if (iammaster) then
      write (*,*) ' Calculating three-center exchange-correlation '
      write (*,*) ' interactions (OLSXC). '
    end if
    nstyles = ideriv
    call gleg (ctheta, ctheta_weights, ntheta_max)
    do looper3a = 1, nspec*nspec*nspec
      if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) &
      &        .eq. my_proc) then
        itmp   = looper3a
        itype3 = 1 + int((itmp - 1)/(nspec*nspec))
        itmp   = itmp - (itype3 - 1)*(nspec*nspec)
        itype2 = 1 + int((itmp - 1)/nspec)
        itmp   = itmp - (itype2 - 1)*nspec
        itype1 = itmp

        rcutoff1 = rcutoffa_max(itype1)
        rcutoff2 = rcutoffa_max(itype2)
        rcutoff3 = rcutoffa_max(itype3)

        atom1 = atom(itype1)
        atom2 = atom(itype2)
        atom3 = atom(itype3)

        what1 = what(itype1)
        what2 = what(itype2)
        what3 = what(itype3)

        dbc = rcutoff1 + rcutoff2
        dna = rcutoff3 + max(rcutoff1,rcutoff2)

        index_max = index_max3c(itype1,itype2)
        do index = 1, index_max
          n1(index) = nleft(itype1,itype2,index)
          l1(index) = lleft(itype1,itype2,index)
          m1(index) = mleft(itype1,itype2,index)
          n2(index) = nright(itype1,itype2,index)
          l2(index) = lright(itype1,itype2,index)
          m2(index) = mright(itype1,itype2,index)
        end do

! If the cutoffs for the atoms are too drastically different, then the grid
! size should be increased. Otherwise, the number of non-zero points may
! be too few.
        max_diff = max(abs(rcutoff1 - rcutoff2), &
        &                    abs(rcutoff2 - rcutoff3), &
        &                    abs(rcutoff3 - rcutoff1))
        if (max_diff .gt. 2.0_dp) then
          write (*,*) ' ************ WARNING ************* '
          write (*,*) ' You have at least two species which have '
          write (*,*) ' rcutoff''s which differ by more than 2.0 '
          write (*,*) ' Angstroms. It is highly advisable that you '
          write (*,*) ' increase the number of mesh points, so as to '
          write (*,*) ' avoid a case where you may end up with many '
          write (*,*) ' zeros, and too few non-zero elements in your '
          write (*,*) ' grid. '
        end if
! xc3c_SN now we call threecenter() with new option 'interaction=3'
        interaction = 3
        ispher = .true.

! Even in case of harris option we need all interactions
        isorpmin = 1
        isorpmax = nssh(itype3)

! Do not parallelize over ispnum, because it involves little work
        ispnum = isorpmax - isorpmin + 1
        call threecenter (itype1, itype2, itype3, index_max, iexc, &
        &                       interaction, nzx, nssh, n1, l1, m1, &
        &                       n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
        &                       atom1, atom2, atom3, what1, what2, what3, &
        &                       dbc, dna, signature, ctheta,  &
        &                       ctheta_weights, isorpmin,  &
        &                       isorpmax, ispnum, iammaster, ispher)

      end if ! end MPI which node
    end do
  end if
! xc3c_SN: end of the added part


! MPI CLEAN UP
  if (iammaster) then
    call pp_end()
    call pot_end()
    call wf_end()
    call xc_end()
  end if
  call Finalize_MPI

! Format Statements
! ======================================================================
100 format (70('='))
102 format (a70)
103 format (2x, ' Hello ', a70,/,' Good to see you! ')
104 format (2x, a70)
106 format (a45)
107 format (' (', 6(i2,',',1x), i2,')')
108 format (' (', 2(i2,',',1x), i2,')')
109 format (' (', 19(i2,',',1x), i2, ') ')
301 format (2x, i2, 9x, ' - Information for this species ')
302 format (2x, a2, 9x, ' - Element ')
303 format (2x, i3, 8x, ' - Nuclear Z ')
304 format (2x, f7.3, 4x, ' - Atomic Mass ')
305 format (2x, i2, 9x, ' - Number of shells; L for each shell ')
306 format (2x, 8(2x,i1))
307 format (2x, i2, 9x, ' - Number of shells; L for each shell ', &
  &                      ' (Pseudopotential) ')
308 format (2x,8(2x,i1))
309 format (8(2x,f5.2), ' - Occupation numbers ')
310 format (8(2x,f5.2), ' - Radial cutoffs ')
311 format (9(2x,a25))
312 format (9(2x,a25))
313 format (2x, f12.5, 2x, ' - Atomic energy ')
314 format ((2x,f5.2),' - Radial cutoffs PP ')

500 format (14x, ' < N L M |  | n l m> ')
501 format (' index = ', i3, 3x, '<',3i2,' |  |',3i2,'>')

end
