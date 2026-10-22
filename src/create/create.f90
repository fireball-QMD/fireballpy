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
use, intrinsic :: iso_fortran_env, only: dp => real64
use :: utils, only: utils_log, UTILS_LOG_INFO, UTILS_LOG_ERROR
use :: xc, only: xc_init, xc_end
use :: wavefunctions, only: wf_init, wf_end
use :: potentials, only: pot_init, pot_end
use :: pseudopotentials, only: pp_init, pp_end
use :: onecenter, only: onecenter_calc, &
  & ONECENTER_XC, &
  & ONECENTER_GOVERLAP1, &
  & ONECENTER_GOVERLAP2
use :: twocenter, only: twocenter_calc, &
  & TWOCENTER_DENS_ATOM, &
  & TWOCENTER_DENS_ONTOP, &
  & TWOCENTER_OVERLAP, &
  & TWOCENTER_VNEUTRAL_ATOM, &
  & TWOCENTER_VNN_ATOM, &
  & TWOCENTER_VNEUTRAL_ONTOP, &
  & TWOCENTER_VNN_ONTOP, &
  & TWOCENTER_VPP, &
  & TWOCENTER_VXC, &
  & TWOCENTER_DIP_Z, &
  & TWOCENTER_DIP_X, &
  & TWOCENTER_DIP_Y, &
  & TWOCENTER_COULOMB
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

! Local Variable Declaration and Description
! ======================================================================
logical :: ispher
integer :: ideriv, iderivmin, iderivmax, iderivtmp, index, index_max, index_maxPP, index_maxsp, &
	   interaction, ispec, isorp, isorpmax, isorpmin, ispnum, istyle, issh, itmp, itype1, itype2, &
	   itype3, nstyles, looper2, looper23, looper24, looper25, looper3, looper3a, minderiv2c, &
	   maxderiv2c, nspec, nssh2, nssh2PP, nzx1, nzx2, itype, lalpha, malpha, nalpha
real(kind=dp) :: fraction, rcutoff1, rcutoff2, rcutoff3, dbc, dna
character(2) :: atom1, atom2, atom3
character(70) :: signature, what1, what2, what3
integer :: index_max2c(nspec_max, nspec_max), index_max2cPP(nspec_max, nspec_max),  &
	   index_max3c(nspec_max, nspec_max), lleft(nspec_max, nspec_max, inter_max), &
	   lleftPP(nspec_max, nspec_max, inter_max), lright(nspec_max, nspec_max, inter_max), &
	   lrightPP(nspec_max, nspec_max, inter_max), mleft(nspec_max, nspec_max, inter_max), &
	   mleftPP(nspec_max, nspec_max, inter_max), mright(nspec_max, nspec_max, inter_max), &
	   mrightPP(nspec_max, nspec_max, inter_max), nleft(nspec_max, nspec_max, inter_max), &
	   nleftPP(nspec_max, nspec_max, inter_max), nright(nspec_max, nspec_max, inter_max), &
	   nrightPP(nspec_max, nspec_max, inter_max), l1(inter_max), l1PP(inter_max), l2(inter_max), &
	   l2PP(inter_max), m1(inter_max), m1PP(inter_max), m2(inter_max), m2PP(inter_max), &
	   n1(inter_max), n1PP(inter_max), n2(inter_max), n2PP(inter_max), n1sp(1), l1sp(1), &
	   m1sp(1), n2sp(1), l2sp(1), m2sp(1), lssh(nspec_max, nsh_max), lsshPP(nspec_max, nsh_max), &
	   nssh(nspec_max), nsshPP(nspec_max, nzx(nspec_max))
real(kind=dp) :: ctheta(ntheta_max), ctheta_weights(ntheta_max), etotatom(nspec_max), &
		 rcutoff(nspec_max, nsh_max), rcutoffa(nspec_max, nsh_max), rcutoffa_max(nspec_max), xmass(nspec_max)
character(2) atom(nspec_max)
character(70) what(nspec_max)

! MPI
integer :: my_proc, nproc
logical :: iammaster, iammpi

! Procedure
! ======================================================================
! Setup the clm coefficients
call setup_clm()

! Initialize MPI
call init_MPI(iammaster, iammpi, my_proc, nproc)

if (iammaster) then
  call utils_log(UTILS_LOG_INFO, "Running basis creation on "//trim(cnproc)//" nodes")
  signature = ' fireballpy '
end if ! end master

!  MPI Broadcast signature
call signature_MPI (signature)

! **********************************************************************
do itype1 = 1, nspec
  do itype2 = 1, nspec
    call mk_index(itype1, itype2, nspec_max, nsh_max, inter_max, &
                  nssh, lssh, nleft, lleft, mleft, nright, &
                  lright, mright, index_max2c, index_max3c)

    ! TODO: this should not happen
    if (index_max3c(itype1,itype2) .gt. inter_max) then
      call utils_log(UTILS_LOG_ERROR, "inter_max too small")
      error stop 1
    end if

    ! For the pseudopotential call a different mk_index routine
    call mk_indexPP(itype1, itype2, nspec_max, nsh_max, &
                    inter_max, nssh, nsshPP, lssh, lsshPP, &
                    nleftPP, lleftPP, mleftPP, nrightPP, &
                    lrightPP, mrightPP, index_max2cPP)
  end do
end do

! Done with setup, now get to work
! ======================================================================
if (iammaster) then
  call wf_init()
  call pot_init()
  call pp_init()

  call utils_log(UTILS_LOG_INFO, "Computing one center integrals")
  ! TODO: mpi
  do itype1 = 1, nspec
    call onecenter_calc(ONECENTER_XC + &
                       ONECENTER_GOVERLAP1 + &
                       ONECENTER_GOVERLAP2, ispec=itype1)
  end do

  call utils_log(UTILS_LOG_INFO, "Computing two center integrals")
  ! TODO: mpi
  do itype1 = 1, nspec
    do itype2 = 1, nspec
      call twocenter_calc(TWOCENTER_DENS_ATOM + &
                         TWOCENTER_DENS_ONTOP + &
                         TWOCENTER_OVERLAP + &
                         TWOCENTER_VNEUTRAL_ATOM + &
                         TWOCENTER_VNN_ATOM + &
                         TWOCENTER_VNEUTRAL_ONTOP + &
                         TWOCENTER_VNN_ONTOP + &
                         TWOCENTER_VPP + &
                         TWOCENTER_VXC + &
                         TWOCENTER_DIP_Z + &
                         TWOCENTER_DIP_X + &
                         TWOCENTER_DIP_Y + &
                         TWOCENTER_COULOMB, ispec=itype1, jspec=itype2)
      call twocenter_calc(TWOCENTER_DENS_ATOM + &
                         TWOCENTER_DENS_ONTOP + &
                         TWOCENTER_OVERLAP, ispec=itype1, jspec=itype2, is_spherical=.true.)
    end do
  end do

end if ! end master

! ======================================================================
! I. Perform three-center calculations
! ======================================================================
! **********************************************************************
!
!  =====>         1b. Three center exchange-correlation matrix element
!                                   (three-center SNXC)
!
! **********************************************************************

ideriv = 0

if (iammaster) then
  call utils_log(UTILS_LOG_INFO, "Computing three center integrals")
  call utils_log(UTILS_LOG_INFO, "  SNXC")
end if
nstyles = ideriv
call gleg (ctheta, ctheta_weights, ntheta_max)
do looper3a = 1, nspec*nspec*nspec
  if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) &
    & .eq. my_proc) then
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

    ! xc3c_SN now we call threecenter() with new option 'interaction=3'
    interaction = 3
    ispher = .false.

    ! Even in case of harris option we need all interactions
    isorpmin = 1
    isorpmax = nssh(itype3)

    ! Do not parallelize over ispnum, because it involves little work
    ispnum = isorpmax - isorpmin + 1
    call threecenter(itype1, itype2, itype3, index_max, &
                     interaction, nzx, nssh, n1, l1, m1, &
                     n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
                     atom1, atom2, atom3, what1, what2, what3, &
                     dbc, dna, signature, ctheta,  &
                     ctheta_weights, isorpmin, &
                     isorpmax, ispnum, iammaster, ispher)

  end if ! end MPI which node
end do
! xc3c_SN: end of the added part


! **********************************************************************
!
!  =====>         2. Three center neutral atom matrix element
!
! **********************************************************************
if (iammaster) then
  call utils_log(UTILS_LOG_INFO, "  charged atom interactions")
end if
call gleg (ctheta, ctheta_weights, ntheta_max)
do looper3a = 1, nspec*nspec*nspec
  if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) .eq. my_proc) then
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

    interaction = 1
    ispher = .false.
    isorpmin = 0
    isorpmax = nssh(itype3)

    ! Do not parallelize over ispnum, because it involves little work
    ispnum = isorpmax - isorpmin + 1
    call threecenter (itype1, itype2, itype3, index_max, &
                     interaction, nzx, nssh, n1, l1, m1, &
                     n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
                     atom1, atom2, atom3, what1, what2, what3, &
                     dbc, dna, signature, &
                     ctheta, ctheta_weights, isorpmin, &
                     isorpmax, ispnum, iammaster, ispher)
  end if ! end MPI which node
end do

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
    call kinetic (itype1, itype2, atom1, atom2, what1, what2, &
                 nzx1, nzx2, rcutoff1, rcutoff2, nssh, lssh, &
                 index_max, n1, l1, m1, n2, l2, m2, signature, &
                 iammaster)
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
                  nssh, lssh, nleft, lleft, mleft, nright, &
                  lright, mright, index_max2c, index_max3c)
  end do
end do

! **********************************************************************
!
!  =====>         1b. Three center exchange-correlation matrix element
!                                   (three-center OLSXC)
!                                  SPHERICAL APPROXIMATION
! **********************************************************************

ideriv = 0

if (iammaster) then
  call utils_log(UTILS_LOG_INFO, "  OLSXC")
end if
nstyles = ideriv
call gleg (ctheta, ctheta_weights, ntheta_max)
do looper3a = 1, nspec*nspec*nspec
  if (mod(looper3a + nspec*nspec*nspec*nstyles, nproc) &
    & .eq. my_proc) then
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

    ! xc3c_SN now we call threecenter() with new option 'interaction=3'
    interaction = 3
    ispher = .true.

    ! Even in case of harris option we need all interactions
    isorpmin = 1
    isorpmax = nssh(itype3)

    ! Do not parallelize over ispnum, because it involves little work
    ispnum = isorpmax - isorpmin + 1
    call threecenter (itype1, itype2, itype3, index_max, &
                     interaction, nzx, nssh, n1, l1, m1, &
                     n2, l2, m2, rcutoff1, rcutoff2, rcutoff3, &
                     atom1, atom2, atom3, what1, what2, what3, &
                     dbc, dna, signature, ctheta,  &
                     ctheta_weights, isorpmin,  &
                     isorpmax, ispnum, iammaster, ispher)

  end if ! end MPI which node
end do
! xc3c_SN: end of the added part


if (iammaster) then
  call pp_end()
  call pot_end()
  call wf_end()
  call xc_end()
end if

! MPI CLEAN UP
call Finalize_MPI

end program create
