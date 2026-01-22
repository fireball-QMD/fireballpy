! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! awp with this program.  If not, see <http://www.gnu.org/licenses/>.


! onecenter.f90
! Program Description
! ==============================================================================
!       This module calculates the one-center integrals for the exchange-
! correlation interactions.
! ==============================================================================
! Original code from Juergen Fritsch
!
! Code rewritten by:
! C. Roldan Pinero
! ==============================================================================
module onecenter
  use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
  use precision, only: wp
  use constants, only: abohr, hartree
  use math, only: math_lstsq
  implicit none
  private
  public :: onecenter_init
  public :: onecenter_calc

  integer :: iexc_, nsh_max_, nspec_, nspec_max_, wfmax_points_
  integer, allocatable :: iderorb_(:), nsshxc_(:), nzx_(:), lsshxc_(:,:)
  real(kind=wp) :: fraction_
  real(kind=wp), allocatable :: dqorb_(:), drr_rho_(:), rcutoffa_max_(:),      &
  & xnocc_(:,:)

contains

  ! TODO: all these things should be read from a module
  integer function onecenter_init(nspec, nspec_max, nsh_max, wfmax_points,     &
  & iexc, fraction, nsshxc, lsshxc, rcutoffa_max, xnocc, dqorb, iderorb,       &
  & drr_rho, nzx)
    integer, intent (in) :: iexc, nsh_max, nspec, nspec_max, wfmax_points,     &
    & iderorb(:), nsshxc(:), nzx(:), lsshxc(:,:)
    real(kind=wp), intent (in) :: fraction, dqorb(:), drr_rho(:),              &
    & rcutoffa_max(:), xnocc(:,:)

    allocate(iderorb_(nspec_max), nsshxc_(nspec_max), nzx_(nspec_max),         &
    & lsshxc_(nspec_max, nsh_max), dqorb_(nspec_max), drr_rho_(nspec_max),     &
    & rcutoffa_max_(nspec_max), xnocc_(nsh_max, nspec_max), stat=onecenter_init)
    if (onecenter_init /= 0) then
      write(stderr, '(a)') '[ERROR] onecenter.f90: failed allocation'
      return
    end if

    iexc_ = iexc
    nsh_max_ = nsh_max
    nspec_ = nspec
    nspec_max_ = nspec_max
    wfmax_points_ = wfmax_points
    fraction_ = fraction
    iderorb_ = iderorb
    nsshxc_ = nsshxc
    nzx_ = nzx
    lsshxc_ = lsshxc
    dqorb_ = dqorb
    drr_rho_ = drr_rho
    rcutoffa_max_ = rcutoffa_max
    xnocc_ = xnocc

    onecenter_init = 0
    return
  end function onecenter_init

  integer function onecenter_calc()
    integer :: ispec
    do ispec = 1, nspec_
      onecenter_calc = onecenter_spec_calc(ispec)
      if (onecenter_calc /= 0) then
        call onecenter_free()
        write(stderr, '(a,i4)') '[ERROR] onecenter.f90: failed calc for ispec =', ispec
        return
      end if
    end do
    call onecenter_free()
    onecenter_calc = 0
    return
  end function onecenter_calc

  integer function onecenter_spec_calc(ispec)
    integer, intent(in) :: ispec
    logical :: first_orb(0:2)
    integer :: irho, issh, jssh, kssh, ilssh, jlssh, iorb, jorb,               &
    & nssh, nssh1, ix, ndq, nnz, npts, nnrho, nrhs, norbs, im
    real(kind=wp) :: dnuxc, dnuxcs, drho, exc, dexcc, factor, rcutoff, rho,    &
    & rhomax, rhomin, rh, rhp, rhpp, vxc, ddq, tmp, xnocc_in(nsh_max_),        &
    & rho1c(wfmax_points_), rhop1c(wfmax_points_), rhopp1c(wfmax_points_)
    character(len=2) :: auxz
    character(len=256) :: errmsg
    logical, allocatable :: iszero(:)
    integer, allocatable :: orb2ssh(:), morbs(:)
    real(kind=wp), allocatable :: xmatt(:,:), eexc(:,:), vvxc(:,:),            &
    & tmpmat(:,:), dq(:,:), exc1crho(:,:,:), nuxc1crho(:,:,:)
    real(kind=wp), external :: psiofr

    nssh = nsshxc_(ispec)
    nssh1 = nssh + 1

    ! Needed for charge corrections
    drho = drr_rho_(ispec)
    rcutoff = rcutoffa_max_(ispec)
    rhomin = 0.0_wp
    rhomax = rcutoff
    nnrho = nint((rhomax - rhomin)/drho) + 1

    ! TODO: this should come from somewhere else
    ! Prepare the increments
    first_orb = .true.
    ndq = 3
    npts = ndq**nssh
    nnz = (nssh*(nssh + 1))/2
    allocate(dq(nssh, ndq))
    do issh = 1, nssh
      ilssh = lsshxc_(ispec, issh)
      if (ilssh == 0) then
        if (first_orb(ilssh)) then
          dq(issh, 1) = -0.10_wp
          dq(issh, 2) = 0.00_wp
          dq(issh, 3) = 0.10_wp
          first_orb(ilssh) = .false.
        else
          dq(issh, 1) = 0.00_wp
          dq(issh, 2) = 0.05_wp
          dq(issh, 3) = 0.10_wp
        end if
      else if (ilssh == 1) then
        if (first_orb(ilssh)) then
          dq(issh, 1) = -0.20_wp
          dq(issh, 2) = 0.00_wp
          dq(issh, 3) = 0.20_wp
          first_orb(ilssh) = .false.
        else
          dq(issh, 1) = 0.00_wp
          dq(issh, 2) = 0.10_wp
          dq(issh, 3) = 0.20_wp
        end if
      else
        if (first_orb(ilssh)) then
          dq(issh, 1) = -0.50_wp
          dq(issh, 2) = 0.00_wp
          dq(issh, 3) = 0.50_wp
          first_orb(ilssh) = .false.
        else
          dq(issh, 1) = 0.00_wp
          dq(issh, 2) = 0.25_wp
          dq(issh, 3) = 0.50_wp
        end if
      end if
    end do

    ! We get the full grid of values first
    allocate(iszero(nnz))
    iszero = .false.
    kssh = 0
    norbs = 0
    do issh = 1, nssh
      ilssh = lsshxc_(ispec, issh)
      norbs = norbs + (2*ilssh + 1)
      do jssh = issh, nssh
        kssh = kssh + 1
        if (ilssh /= lsshxc_(ispec, jssh)) iszero(kssh) = .true.
      end do
    end do
    allocate(eexc(npts, nnz), vvxc(npts, nnz))
    eexc = 0.0_wp
    vvxc = 0.0_wp
    allocate(xmatt(nssh1, npts))
    do ix = 1, npts
      ! Set charges
      xmatt(1, ix) = 1.0_wp
      do issh = 1, nssh
        ddq = dq(issh, 1 + mod(ix - 1, ndq**issh)/ndq**(issh - 1))
        xnocc_in(issh) = xnocc_(issh, ispec) + ddq
        xmatt(issh + 1, ix) = ddq
      end do

      ! Obtain the density and respective derivatives needed for evaluating the
      ! exchange-correlation interactions (LDA or GGA).
      ! We have to avoid change xnocc_in !!
      rho1c = 0.0_wp
      rhop1c = 0.0_wp
      rhopp1c = 0.0_wp
      call rho1c_store (ispec, nsh_max_, nssh, 0.0_wp, 1, drho, rcutoff,       &
      &                 xnocc_in, 1, wfmax_points_, rho1c, rhop1c, rhopp1c)

      ! Integrals <i|exc(i)|i> and <i.nu|mu(i)|i.nu'>
      do irho = 1, nnrho
        rho = rhomin + real(irho - 1, kind=wp)*drho
        factor = 0.66666666666666666667_wp*drho
        if (mod(irho, 2) == 0) factor = 2.0_wp*factor
        if (irho == 1 .or. irho == nnrho) factor = 0.5_wp*factor

        ! Compute the exchange correlation potential
        rho = rho/abohr
        rh = rho1c(irho)*abohr**3
        rhp = rhop1c(irho)*abohr**4
        rhpp = rhopp1c(irho)*abohr**5
        call get_potxc1c (iexc_, fraction_, rho, rh, rhp, rhpp, exc, vxc,      &
        &                 dnuxc, dnuxcs, dexcc)
        vxc = hartree*vxc
        exc = hartree*exc
        rho = rho*abohr

        kssh = 0
        do issh = 1, nssh
          do jssh = issh, nssh
            kssh = kssh + 1
            if (iszero(kssh)) cycle
            tmp = psiofr(ispec, issh, rho)*psiofr(ispec, jssh, rho)*factor*rho*rho
            eexc(ix, kssh) = eexc(ix, kssh) + tmp*exc
            vvxc(ix, kssh) = vvxc(ix, kssh) + tmp*vxc
          end do
        end do
      end do
    end do
    deallocate(dq)

    ! Perform the fits
    nrhs = 2
    allocate(tmpmat(npts, nrhs*nnz))
    tmpmat(:, 1:nnz) = eexc
    tmpmat(:, nnz+1:2*nnz) = vvxc
    deallocate(eexc, vvxc)
    onecenter_spec_calc = math_lstsq(xmatt, tmpmat, is_a_trans=.true.)
    deallocate(xmatt)
    if (onecenter_spec_calc /= 0) then
      deallocate(tmpmat, iszero)
      write (stderr, '(a,i4)') '[ERROR] onecenter.f90: call to `(x)gelst` failed with info =', &
      &                        onecenter_spec_calc
      return
    end if

    ! Prepare for output
    allocate(exc1crho(0:nssh,nssh,nssh), nuxc1crho(0:nssh,nssh,nssh))
    exc1crho = 0.0_wp
    nuxc1crho = 0.0_wp
    do ix = 0, nssh
      kssh = 0
      do issh = 1, nssh
        do jssh = issh, nssh
          kssh = kssh + 1
          if (iszero(kssh)) cycle
          exc1crho(ix, issh, jssh) = tmpmat(ix+1, kssh)
          exc1crho(ix, jssh, issh) = tmpmat(ix+1, kssh)
          nuxc1crho(ix, issh, jssh) = tmpmat(ix+1, kssh + nnz)
          nuxc1crho(ix, jssh, issh) = tmpmat(ix+1, kssh + nnz)
        end do
      end do
    end do
    deallocate(tmpmat, iszero)

    ! Useful mapping
    allocate(orb2ssh(norbs), morbs(norbs))
    iorb = 0
    do issh = 1, nssh
      im = -lsshxc_(ispec, issh)
      do ilssh = 1, (2*lsshxc_(ispec, issh) + 1)
        iorb = iorb + 1
        orb2ssh(iorb) = issh
        morbs(iorb) = im
        im = im + 1
      end do
    end do

    ! Write log
    write (auxz,'(i2.2)') nzx_(ispec)
    write (stdout,'(a)') 'Writing output to: coutput/onecenter_xc.'//auxz//'.dat'

    ! Write output
    open (unit=360, file='coutput/onecenter_xc.'//auxz//'.dat', status='unknown')
    write (360,'(2i4)') nssh, norbs
    do ix = 0, nssh
      do iorb = 1, norbs
        do jorb = 1, norbs
          if (morbs(iorb) == morbs(jorb)) then
            write (360,'(ES15.7)', advance='no') nuxc1crho(ix, orb2ssh(iorb), orb2ssh(jorb))
          else
            write (360,'(ES15.7)', advance='no') 0.0_wp
          end if
        end do
        write (360,'(a)') ''
      end do
      write (360,'(a)') ''
    end do
    do ix = 0, nssh
      write (360,'(100ES15.7)') (exc1crho(ix, issh, issh), issh = 1, nssh)
      write (360,'(a)') ''
    end do
    close(360)
    deallocate(orb2ssh, morbs, exc1crho, nuxc1crho)

    onecenter_spec_calc = 0
    return
  end function onecenter_spec_calc

  subroutine onecenter_free()
    deallocate(iderorb_, nsshxc_, nzx_, lsshxc_, dqorb_, drr_rho_, rcutoffa_max_, xnocc_)
  end subroutine onecenter_free
end module onecenter
