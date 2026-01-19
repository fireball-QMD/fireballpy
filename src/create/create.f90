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
! geometry. This new creator includes the following subroutines.
!
! coulomb.f -- coulomb integrals for the hartree energy
! kinetic.f -- kinetic energy matrix elements (pure 2-center)

! These two cases have been replaced by the Kleinman-Bylander
! pseupotential so that everything is now only two center.
! nl_atom.f --neutral atom, nonlocal pseudopotential degenerate 3-center
! nl_ontop.f --neutral atom, on top case degenerate 3-center

! These six cases are now combined into one general two-center
! integration routine called twocenter.f
!  dipole.f -- dipole interactions
!  overlap.f -- overlap between 2 centers
!  na_ontop.f -- neutral atom on top degenerate 3 center case
!  na_atom.f -- neutral atom "atom" case deg. 3 center
!  xc_atom.f -- exchange correlation -- atom case  deg 3 center
!  xc_ontop.f -- exchange correlation -- ontop case deg 3 center

!  xcna.f XC -- true three center case
!  bcna.f BC-NA -- true three center case
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

! ======================================================================
!
! Program Declaration
! ======================================================================
program create
  use precision, only: wp
  use onecenter, only: onecenter_init, onecenter_calc
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
  real(kind=wp) abohr
  parameter (abohr = 0.529177249d0)

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
  real(kind=wp) ctheta (ntheta_max)
  real(kind=wp) ctheta_weights (ntheta_max)

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

  real(kind=wp) fraction
  real(kind=wp) max_diff
  real(kind=wp) rcutoff1
  real(kind=wp) rcutoff2
  real(kind=wp) rcutoff3

! distances needed for three center integrals
  real(kind=wp) dbc                      ! maximal bond charge distance
  real(kind=wp) dna                      ! maximal neutral atom distance

  real(kind=wp) etotatom (nspec_max)
  real(kind=wp) rcutoff (nspec_max, nsh_max) ! cutoff radius in bohr
  real(kind=wp) rcutoffa (nspec_max, nsh_max)! cutoff radius in angstroms
  real(kind=wp) rcutoffa_max (nspec_max)     ! cutoff radius in angstroms
  real(kind=wp) xmass (nspec_max)

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
    write (*,*) '          Latest version Jan 10, 2005 '
    write (*,*) '          See Copyright information: '
    write (*,*) '            !!!Proprietory Code!!! '
    write (*,*) '       Usable only with permission from '
    write (*,*) '      the Fireball executive committee. '
    write (*,*) '  This program is NOT, under any circumstances '
    write (*,*) '    to be transfered to an unauthorized user. '
    write (*,100)
    write (*,*) '  '

! Identification of the user:
!         write (*,*) '  '
!         write (*,*) ' Please insert your name and other messages. '
!         write (*,*) '  '
    if (.not. iammpi) then
!          read (*,102) signature
      signature = ' fireballpy '
    else
! You might want to edit this to put your own name in here
      signature = ' MPI superhero '
    endif
    write (*,*) '  '
    write (*,103) signature
    write (*,*) '  '
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
    if (iammaster) then
      write (*,*) ' The exchange-correlation option that you chose '
      write (*,*) ' has not been implemented into creator yet. '
      write (*,*) ' Choose a different one, and restart. '
    end if
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
    write (*,*) '  '
    write (*,*) '  '
    write (*,*) ' Now we are writing the following to the '
    write (*,*) ' info.dat file. '
    write (*,*) '  '

! First write to the screen
    do ispec = 1, nspec
      write (*,100)
      write (*,301) ispec
      write (*,302) atom(ispec)
      write (*,303) nzx(ispec)
      write (*,304) xmass(ispec)
      write (*,305) nssh(ispec)
      write (*,306) (lssh(ispec,issh), issh = 1, nssh(ispec))
      write (*,307) nsshPP(ispec)
      write (*,308) (lsshPP(ispec,issh), issh = 1, nsshPP(ispec))
!jel-PP
      write (*,314) (rcPP(ispec))

      write (*,309) (xnocc(issh,ispec), issh = 1, nssh(ispec))
      write (*,310) (rcutoff(ispec,issh), issh = 1, nssh(ispec))
      write (*,311) (wavefxn(ispec,issh), issh = 1, nssh(ispec))
      write (*,312) (napot(ispec,issh), issh = 0, nssh(ispec))
      write (*,313) etotatom(ispec)

      write (*,100)
    end do

! Next write to the info.dat file.
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
    write (*,*) '  '

! Done with setup, now get to work
! ======================================================================

! O. Compute exchange one-center case.
! **********************************************************************
!
!  =====>         4.5  Dipole???
!
! **********************************************************************
    if (V_intra_dip .eq. 1) then
      do itype = 1,nspec
        write(*,*) 'itype= ', itype !Ankais
        call onecentervdip (nsh_max, nspec, nspec_max, itype, &
        &               nssh, lssh, drr_rho, rcutoffa_max, what, signature)
      end do !end do itype
    end if !end if (V_intra_dip .eq. 1) then
! ======================================================================
! Only do this on the master
    if (imuxc1c .eq. 1) then
! JOM-add
      call goverlap1c (nspec, nspec_max, nsh_max, wfmax_points,  &
      &                      nssh, lssh, rcutoffa_max, what,  &
      &                      signature, drr_rho)
! JOM-end
      if (onecenter_init(nspec, nspec_max, nsh_max, wfmax_points, iexc,        &
      &                  fraction, nsshxc, lsshxc, rcutoffa_max, xnocc, dqorb, &
      &                  iderorb, drr_rho, nzx) /= 0) stop 1
      if (onecenter_calc() /= 0) stop 1
    end if

! JPL 1999 Exact exchange interactions.
    if (imuxc1c .eq. 1 .and. iexc .eq. 12) then
      call x_1c (nsh_max, nspec, nspec_max, fraction, nssh, lssh, &
      &               drr_rho, rcutoffa_max, what, signature)
    end if
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
!  =====>         1. Three center exchange-correlation matrix element
!
! **********************************************************************

        ideriv = 0

        if (ibcxc .eq. 1 .and. ixc_opt .eq. 0) then
         if (iammaster) then
          write (*,*) ' Calculating three-center exchange-correlation '
          write (*,*) ' interactions (Horsfield). '
         end if
         if (idogs .eq. 1) then
          iderivmin = 1
          iderivmax = 6
         else
          iderivmin = 0
          iderivmax = 0
         end if
         if (iharris .eq. 1) iderivmin = 0
         ideriv = iderivmax - iderivmin + 1

! If a GGA was chosen, then inform the user that GGA's cannot currently
! be used for the three-center case.  Then by default set iexc to the
! Ceperley-Alder/Perdew-Zunger form of LDA (iexc = 3).
         if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
     &       .or. iexc .eq. 9 .or. iexc .eq. 10) then
          if (iammaster) then
           write (*,*) '  '
           write (*,*) ' You have chosen to perform a GGA type of '
           write (*,*) ' exchange-correlation interaction. However, '
           write (*,*) ' currently this capability does not exist '
           write (*,*) ' for three-center interactions. By default, '
           write (*,*) ' we set the three-center interactions to the '
           write (*,*) ' LDA limit - Ceperley-Alder/Perdew-Zunger. '
          end if
          iexc_new = 3
         else
          iexc_new = iexc
         end if
         nstyles = ideriv

! if ibcxc .eq. 0, then this loop is skipped since nstyles = 0
         call gleg (ctheta, ctheta_weights, ntheta_max)
         do looper3 = 1, nspec*nspec*nspec*nstyles
          if (mod(looper3,nproc) .eq. my_proc) then
           itmp   = looper3
           itype3 = 1 + int((itmp - 1)/(nspec*nspec*nstyles))
           itmp   = itmp - (itype3 - 1)*(nspec*nspec*nstyles)
           itype2 = 1 + int((itmp - 1)/(nspec*nstyles))
           itmp   = itmp - (itype2 - 1)*(nspec*nstyles)
           itype1 = 1 + int((itmp - 1)/nstyles)
           itmp   = itmp - (itype1 - 1)*nstyles
           istyle = itmp

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
           if (max_diff .gt. 2.0d0) then
            write (*,*) ' ************ WARNING ************* '
            write (*,*) ' You have at least two species which have '
            write (*,*) ' rcutoff''s which differ by more than 2.0 '
            write (*,*) ' Angstroms. It is highly advisable that you '
            write (*,*) ' increase the number of mesh points, so as to '
            write (*,*) ' avoid a case where you may end up with many '
            write (*,*) ' zeros, and too few non-zero elements in your '
            write (*,*) ' grid. '
           end if

           interaction = 2
           ispher = .false.
! GGA's cannot currently be used for the three-center case,
! thus use iexc_new instead of iexc.

           iderivtmp = istyle - 1 + iderivmin
           call threecenter (itype1, itype2, itype3, index_max, &
     &                       iexc_new, interaction, nzx, nssh, n1, l1, &
     &                       m1, n2, l2, m2, rcutoff1, rcutoff2, &
     &                       rcutoff3, atom1, atom2, atom3, what1, &
     &                       what2, what3, dbc, dna, signature, &
     &                       ctheta, ctheta_weights, iderivtmp, &
     &                       iderivtmp, 1, iammaster, ispher)
          end if !MPI which node end if
         end do
        end if

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
           if (max_diff .gt. 2.0d0) then
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
           if (max_diff .gt. 2.0d0) then
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


! **********************************************************************
!
!  =====>         2. Overlap
!
! **********************************************************************
          if (iswitch(1) .eq. 1) then
           interaction = 1
           isorp = 0
           ideriv = 0
           ispher = .false.
           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                     nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)

          end if

! **********************************************************************
!
!  =====>         3. Non-local
!
! **********************************************************************
          if (iswitch(4) .eq. 1) then
           index_maxPP = index_max2cPP(itype1,itype2)
           do index = 1, index_maxPP
            n1PP(index) = nleftPP(itype1,itype2,index)
            l1PP(index) = lleftPP(itype1,itype2,index)
            m1PP(index) = mleftPP(itype1,itype2,index)
            n2PP(index) = nrightPP(itype1,itype2,index)
            l2PP(index) = lrightPP(itype1,itype2,index)
            m2PP(index) = mrightPP(itype1,itype2,index)
           end do
           nssh2PP = nsshPP(itype2)

           interaction = 4
           isorp = 0
           ideriv = 0
           ispher = .false.
           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcPP(itype2), nssh2PP, &
     &                     nznl, nrhonl, nddnl, index_maxPP, n1PP, l1PP, &
     &                     m1PP, n2PP, l2PP, m2PP, signature, iammaster, &
     &                     ispher)
          end if

! **********************************************************************
!
!  =====>         4. Exchange-correlation over counting correction
!
! **********************************************************************
! We calculate  (n1+n2)*(exc(1+2)-muxc(1+2)) - n1*(exc(1)-xcmu(1))
!                                            - n2*(exc(2)-xcmu(2))
          if (iswitch(7) .eq. 1) then
           interaction = 7            ! atom/atom
           isorp = 0

! This is not a matrix element, but an over-counting correction to the
! exchange-correlation interaction, so set index_max = 1
           index_maxsp = 1            ! only one term
           n1sp(1) = 1
           l1sp(1) = 0
           m1sp(1) = 0
           n2sp(1) = 1
           l2sp(1) = 0
           m2sp(1) = 0

! Here is the loop over the key.
! Only do derivatives for idogs .eq. 1
           if (idogs .eq. 1) then
            minderiv2c = 1
            maxderiv2c = 4
           else
            minderiv2c = 0
            maxderiv2c = 0
           end if
           if (iharris .eq. 1) minderiv2c = 0
           ispher = .false.
           do ideriv = minderiv2c, maxderiv2c
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, &
     &                      nzxco, nrhoxco, nddxco, index_maxsp, n1sp, &
     &                      l1sp, m1sp, n2sp, l2sp, m2sp, signature, &
     &                      iammaster, ispher)
           end do
          end if


! **********************************************************************
!
!  =====>         5. Dipole
!
! **********************************************************************
          if (iswitch(8) .eq. 1) then
           interaction = 8                         ! z-dipole
           isorp = 0
           ideriv = 0
           ispher = .false.
           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzd, &
     &                     nrhod, nddd, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)

! Worry about these two cases later - they are needed for IR stuff
           interaction = 9                         ! y-dipole
           interaction = 10                        ! x-dipole
          end if

! **********************************************************************
!
!  =====>         6. Short-range Coulomb
!
! **********************************************************************
          if (iswitch(11) .eq. 1) then
           interaction = 11
           isorp = 0
           ideriv = 0
           ispher = .false.
           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzd, &
     &                     nrhod, nddd, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)
          end if

! **********************************************************************
!
!  =====>         7. Extended Hubbard
!                    nu12 = n1*n2 *nu(n1+n2)
!
! **********************************************************************
          if (inuxc2c .eq. 1) then
           interaction = 12
           isorp = 0
           ideriv = 0

! Important note. We have ONLY programmed in ldaxc.f the nu potential
! (nu = dmu/dn) for iexc = 3. So even if you do a different LDA, or even GGA,
! we will use Ceperly/Alder (iexc = 3) for the nu part.
           iexc_new = 3
           ispher = .false.
           call twocenter (interaction, isorp, ideriv, iexc_new, &
     &                     fraction, itype1, itype2, atom1, atom2, &
     &                     what1, what2, nzx1, nzx2, rcutoff1, rcutoff2, &
     &                     nssh2, nzeh, nrhoeh, nddeh, index_max, n1, &
     &                     l1, m1, n2, l2, m2, signature, iammaster, &
     &                     ispher)
          end if

! **********************************************************************
!
!  =====>         8. Extended Hubbard spin dependent interaction
!                    nu12 = n1*n2 *nu(n1+n2), only implemented
!                    for iexc=11
! **********************************************************************
          if (isnuxc2c .eq. 1) then
           if (iexc .eq. 11) then
            interaction = 13
            isorp = 0
            ideriv = 0
            ispher = .false.
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2,  &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzeh, &
     &                      nrhoeh, nddeh, index_max, n1, l1, m1, n2,  &
     &                      l2, m2, signature, iammaster, ispher)
           else
            write(*,*) ' This option is currently implemented only '
            write(*,*) ' for the LSDAVWN exchange correlation model '
           end if
          end if
         end if ! MPI which node end if
        end do

! ======================================================================
! III. Perform degenerate three-center calculations
! ======================================================================
        do looper23 = 1, nspec*nspec
         if (mod(looper23,nproc) .eq. my_proc) then
          itmp   = looper23
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
!  =====>         1. Neutral atom/ontop matrix elements
!
! **********************************************************************
          if (iswitch(2) .eq. 1) then
           interaction = 2            ! atom/ontop
           ispher = .false.
           do ideriv = 1, 2           ! 1 => left, 2 => right

! We do not calculate charge interactions if idogs .ne. 1
            if (idogs .eq. 1) then
             isorpmin = 1
             if (ideriv .eq. 1) isorpmax = nssh(itype1)
             if (ideriv .eq. 2) isorpmax = nssh(itype2)
            else
             isorpmin = 0
             isorpmax = 0
            end if
            if (iharris .eq. 1) isorpmin = 0

            do isorp = isorpmin, isorpmax
             call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                       itype1, itype2, atom1, atom2, what1, what2, &
     &                       nzx1, nzx2, rcutoff1, rcutoff2, nssh2, &
     &                       nznao, nrhonao, nddnao, index_max, n1, &
     &                       l1, m1, n2, l2, m2, signature, iammaster, &
     &                       ispher)
            end do
           end do
          end if

! **********************************************************************
!
!  =====>         2. Average Density  atom/ontop
!
! **********************************************************************
          if (iswitch(0) .eq. 1) then
           ispher = .false.
           interaction = 0
           do ideriv = 1,2

! Even in case of harris option we need all interactions
            isorpmin = 1
            if (ideriv .eq. 1) isorpmax = nssh(itype1)
            if (ideriv .eq. 2) isorpmax = nssh(itype2)

            do isorp = isorpmin, isorpmax
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                      nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                      m2, signature, iammaster, ispher)
              end do
           end do
! McWeda charge transfer correction terms
           interaction = 14
           do ideriv = 1,2

! Even in case of harris option we need all interactions
            isorpmin = 1
            if (ideriv .eq. 1) isorpmax = nssh(itype1)
            if (ideriv .eq. 2) isorpmax = nssh(itype2)

            do isorp = isorpmin, isorpmax
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                      nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                      m2, signature, iammaster, ispher)
            end do
           end do


          end if


! **********************************************************************
!
!  =====>         3. Exchange-correlation atom/ontop matrix elements
!
! **********************************************************************
! We calculate <phi_1|vxc(n1+n2)|phi_2> for the ontop contibution.
! The catch comes in when we compute derivatives. We compute
! neutral,neutral for ideriv1. For other ideriv's we have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 fopr ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
          if (iswitch(5) .eq. 1) then
           interaction = 5            ! atom/ontop

! Only do derivatives for idogs .eq. 1
           if (idogs .eq. 1) then
            minderiv2c = 1
            maxderiv2c = 4
           else
            minderiv2c = 0
            maxderiv2c = 0
           end if
           if (iharris .eq. 1) minderiv2c = 0
           ispher = .false.
           do ideriv = minderiv2c, maxderiv2c
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, &
     &                      nzxco, nrhoxco, nddxco, index_max, n1, &
     &                      l1, m1, n2, l2, m2, signature, iammaster, &
     &                      ispher)
           end do

! JPL 1999 Calculate exact exchange if requested.
           if (iexc .eq. 12) then

! The exchange interaction involves a sum (alpha) over all the orbitals of the
! atom in the potential.  The orbitals are taken as those for the atom on
! the "ket" or left side of the matrix element. Store each term in the sum in
! a separate file.
            isorp = 0
            do nalpha = 1, nssh(itype2)
             lalpha = lssh(itype2,nalpha)
             do malpha = -lalpha, lalpha
              isorp = isorp + 1
              call xontopl_2c (nspec_max, isorp, fraction, &
     &                         nssh, itype1, itype2, atom1, atom2, &
     &                         what1, what2, nzx1, nzx2, rcutoff1, &
     &                         rcutoff2, nzexo, nrhoexo, nddexo, &
     &                         index_max, inter_max, nalpha, lalpha, &
     &                         malpha, n1, l1, m1, n2, l2, m2, &
     &                         signature, iammaster)
              call xontopr_2c (nspec_max, isorp, fraction, &
     &                         nssh, itype1, itype2, atom1, atom2, &
     &                         what1, what2, nzx1, nzx2, rcutoff1, &
     &                         rcutoff2, nzexo, nrhoexo, nddexo, &
     &                         index_max, inter_max, nalpha, lalpha, &
     &                         malpha, n1, l1, m1, n2, l2, m2, &
     &                         signature, iammaster)
             end do
            end do
           end if
          end if

! For the atom/atom interactions, the number of non-zero matrix elements
! is dependent entirely on one atom. Both wavefunctions are located at
! the same site.
          index_max = index_max2c(itype1,itype1)
          do index = 1, index_max
           n1(index) = nleft(itype1,itype1,index)
           l1(index) = lleft(itype1,itype1,index)
           m1(index) = mleft(itype1,itype1,index)
           n2(index) = nright(itype1,itype1,index)
           l2(index) = lright(itype1,itype1,index)
           m2(index) = mright(itype1,itype1,index)
          end do

! **********************************************************************
!
!  =====>         4. Neutral atom/atom matrix elements
!
! **********************************************************************
          if (iswitch(3) .eq. 1) then
           interaction = 3            ! atom/atom
           ideriv = 0
           ispher = .false.
! We do not calculate charge interactions if idogs .ne. 1
           if (idogs .eq. 1) then
            isorpmin = 1
            isorpmax = nssh(itype2)
           else
            isorpmin = 0
            isorpmax = 0
           end if
           if (iharris .eq. 1) isorpmin = 0
           do isorp = isorpmin, isorpmax
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, &
     &                      nznat, nrhonat, nddnat, index_max, n1, &
     &                      l1, m1, n2, l2, m2, signature, iammaster, &
     &                      ispher)
           end do
          end if
! **********************************************************************
!
!  =====>         5. Average Density  atom/atom
!
! **********************************************************************
          if (iswitch(0) .eq. 1) then
           ispher = .false.
           interaction = 0
           ideriv = 0

! Even in case of Harris option we need all interactions
           isorpmin = 1
           isorpmax = nssh(itype2)

           do isorp = isorpmin, isorpmax
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                      nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                      m2, signature, iammaster, ispher)
           end do
          end if

! **********************************************************************
!
!  =====>         6. Exchange-correlation atom/atom matrix elements
!
! **********************************************************************
! We calculate vxc(n1+n2) - vxc(n1) matrix elements.
          if (iswitch(6) .eq. 1) then
           interaction = 6            ! atom/atom
           isorp = 0

! Only do derivatives for idogs .eq. 1
           if (idogs .eq. 1) then
            minderiv2c = 1
            maxderiv2c = 4
           else
            minderiv2c = 0
            maxderiv2c = 0
           end if
           if (iharris .eq. 1) minderiv2c = 0
           ispher = .false.
           do ideriv = minderiv2c, maxderiv2c
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, &
     &                      nzxca, nrhoxca, nddxca, index_max, n1, &
     &                      l1, m1, n2, l2, m2, signature, iammaster, &
     &                      ispher)
           end do

! Calculate exact exchange for doing hybrid exchange interactions.
           if (iexc .eq. 12) then

! The exchange interaction involves a sum (alpha) over all the orbitals of the
! atom in the potential.  The orbitals are taken as those for the atom on
! the "ket" or left side of the matrix element. Store each term in the sum in
! a separate file.
            isorp = 0
            do nalpha = 1, nssh(itype2)
             lalpha = lssh(itype2,nalpha)
             do malpha = -lalpha, lalpha
              isorp = isorp + 1
              call xatom_2c (nspec_max, isorp, fraction, &
     &                       itype1, itype2, atom1, atom2, what1, &
     &                       what2, nzx1, nzx2, rcutoff1, rcutoff2, &
     &                       nzexa, nrhoexa, nddexa, index_max, &
     &                       inter_max, nalpha, lalpha, malpha, n1, l1, &
     &                       m1, n2, l2, m2, signature, iammaster)
             end do
            end do
           end if

! Calculate the exact exchange-type interactions, etc. for the orbital
! occupancy method of Pou et al. PRB 62:4309 (2000).
           if (ioomethod .eq. 1) then
           end if
          end if

         end if ! MPI which node end if
        end do

! ======================================================================
! II. Compute two center cases (Y dipole, X Dipole.
! ======================================================================

        do itype1 = 1, nspec
          do itype2 = 1, nspec
            call mk_indexDipY (itype1, itype2, nspec_max, nsh_max, &
     &              inter_max,nssh, lssh, nleft, lleft, mleft, nright, &
     &              lright, mright, index_max2c, index_max3c)
          enddo
        enddo


        do looper24 = 1, nspec*nspec
         if (mod(looper24,nproc) .eq. my_proc) then
          itmp   = looper24
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


          if (iswitch(9) .eq. 1) then


           interaction = 9                         ! y-dipole
           isorp = 0
           ideriv = 0
           ispher = .false.

           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzd, &
     &                     nrhod, nddd, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)

          end if
         end if ! MPI which node end if
        end do



        do itype1 = 1, nspec
          do itype2 = 1, nspec
            call mk_indexDipX (itype1, itype2, nspec_max, nsh_max, &
     &              inter_max,nssh, lssh, nleft, lleft, mleft, nright, &
     &              lright, mright, index_max2c, index_max3c)
          enddo
        enddo


        do looper25 = 1, nspec*nspec
         if (mod(looper25,nproc) .eq. my_proc) then
          itmp   = looper25
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


          if (iswitch(10) .eq. 1) then


           interaction = 10                        ! x-dipole
           isorp = 0
           ideriv = 0
           ispher = .false.

           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzd, &
     &                     nrhod, nddd, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)

          end if
         end if ! MPI which node end if
        end do



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
           if (max_diff .gt. 2.0d0) then
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


! ======================================================================
! I.Sp. Compute two center cases (overlap, average density)
! with spheric approxomation.
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
!  =====>         2. Average Density (Spherical)
!                        only ON TOP cases
! **********************************************************************
          if (iswitch(0) .eq. 1) then
           ispher = .true.
           interaction = 0
           do ideriv = 1,2

! Even in case of Harris option we need all interactions
            isorpmin = 1
!            if (ideriv .eq. 0) isorpmax = nssh(itype2)
            if (ideriv .eq. 1) isorpmax = nssh(itype1)
            if (ideriv .eq. 2) isorpmax = nssh(itype2)

            do isorp = isorpmin, isorpmax
             call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                      nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                      m2, signature, iammaster, ispher)
            end do
           end do
          end if

! **********************************************************************
!
!  =====>         2. Overlap (Spherical)
!
! **********************************************************************
          if (iswitch(1) .eq. 1) then
           interaction = 1
           isorp = 0
           ideriv = 0
           ispher = .true.
           call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                     itype1, itype2, atom1, atom2, what1, what2, &
     &                     nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                     nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                     m2, signature, iammaster, ispher)

          end if


          index_max = index_max2c(itype1,itype1)
          do index = 1, index_max
           n1(index) = nleft(itype1,itype1,index)
           l1(index) = lleft(itype1,itype1,index)
           m1(index) = mleft(itype1,itype1,index)
           n2(index) = nright(itype1,itype1,index)
           l2(index) = lright(itype1,itype1,index)
           m2(index) = mright(itype1,itype1,index)
          end do

! **********************************************************************
!
!  =====>         5. Average Density  (Spherical)
!                           Atom_case
!
! **********************************************************************
          if (iswitch(0) .eq. 1) then
           ispher = .true.
           interaction = 0
           ideriv = 0

! Even in case of Harris option we need all interactions
           isorpmin = 1
           isorpmax = nssh(itype2)

           do isorp = isorpmin, isorpmax
            call twocenter (interaction, isorp, ideriv, iexc, fraction, &
     &                      itype1, itype2, atom1, atom2, what1, what2, &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh2, nzs, &
     &                      nrhos, ndds, index_max, n1, l1, m1, n2, l2, &
     &                      m2, signature, iammaster, ispher)
           end do
          end if



         end if ! MPI which node end if
        end do

! MPI CLEAN UP
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
