! TODO: E en shells, V tienen que ser cero los ortogonales con misma L

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


! onecenterxc.f90
! Program Description
! ==============================================================================
!       This routine calculates the one-center integrals for the exchange-
! correlation interactions.
! ==============================================================================
! Original code from Juergen Fritsch
!
! Code rewritten by:
! C. Roldan Pinero
! ==============================================================================
subroutine onecenterxc (nspec, nspec_max, nsh_max, wfmax_points,               &
&                       iexc, fraction, nsshxc, lsshxc,                        &
&                       rcutoffa_max, xnocc, dqorb, iderorb,                   &
&                       drr_rho, nzx)
  use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
  use precision, only: wp
  use constants, only: abohr, hartree
  use math, only: math_lstsq
  implicit none

  ! Argument Declaration and Description
  ! ============================================================================
  integer, intent (in) :: iexc, nsh_max, nspec, nspec_max, wfmax_points,       &
  & iderorb(nspec_max), nsshxc(nspec_max), nzx(nspec_max),                     &
  & lsshxc(nspec_max, nsh_max)
  real(kind=wp), intent (in) :: fraction, dqorb(nspec_max),                    &
  & drr_rho(nspec_max), rcutoffa_max(nspec_max), xnocc(nsh_max, nspec_max)

  ! Local Variable Declaration and Description
  ! ============================================================================
  integer :: irho, issh, jssh, kssh, in1, ilssh, jlssh, iorb, jorb,            &
  & nssh, nssh1, ix, ndq, nnz, npts, nnrho, nrhs, info, norbs, im
  real(kind=wp) :: dnuxc, dnuxcs, drho, exc, dexcc, factor, rcutoff, rho,      &
  & rhomax, rhomin, rh, rhp, rhpp, vxc, ddq, tmp, xnocc_in(nsh_max),           &
  & rho1c(wfmax_points), rhop1c(wfmax_points), rhopp1c(wfmax_points)
  character(len=2) :: auxz
  character(len=256) :: errmsg
  logical, allocatable :: iszero(:)
  integer, allocatable :: orb2ssh(:), morbs(:)
  real(kind=wp), allocatable :: xmatt(:,:), eexc(:,:), vvxc(:,:), tmpmat(:,:), &
  & dq(:,:), exc1crho(:,:,:), nuxc1crho(:,:,:)
  real(kind=wp), external :: psiofr

  ! Procedure
  ! ============================================================================
  do in1 = 1, nspec
    nssh = nsshxc(in1)
    nssh1 = nssh + 1

    ! Needed for charge corrections
    drho = drr_rho(in1)
    rcutoff = rcutoffa_max(in1)
    rhomin = 0.0_wp
    rhomax = rcutoff
    nnrho = nint((rhomax - rhomin)/drho) + 1

    ! Prepare the increments
    ndq = 3
    npts = ndq**nssh
    nnz = (nssh*(nssh + 1))/2
    allocate(dq(nssh, ndq))
    do issh = 1,nssh
      ilssh = lsshxc(in1, issh)
      if (ilssh == 0) then
        if (issh == 1) then
          dq(issh,1) = -0.10_wp
          dq(issh,2) = 0.00_wp
          dq(issh,3) = 0.10_wp
        else
          dq(issh,1) = 0.00_wp
          dq(issh,2) = 0.05_wp
          dq(issh,3) = 0.10_wp
        end if
      else if (ilssh == 1) then
        if (issh == 2) then
          dq(issh,1) = -0.20_wp
          dq(issh,2) = 0.00_wp
          dq(issh,3) = 0.20_wp
        else
          dq(issh,1) = 0.00_wp
          dq(issh,2) = 0.10_wp
          dq(issh,3) = 0.20_wp
        end if
      else
        if (issh == 3) then
          dq(issh,1) = -0.50_wp
          dq(issh,2) = 0.00_wp
          dq(issh,3) = 0.50_wp
        else
          dq(issh,1) = 0.00_wp
          dq(issh,2) = 0.25_wp
          dq(issh,3) = 0.50_wp
        end if
      end if
    end do

    ! We get the full grid of values first
    allocate(iszero(nnz))
    iszero = .false.
    kssh = 0
    norbs = 0
    do issh = 1, nssh
      ilssh = lsshxc(in1, issh)
      norbs = norbs + (2*ilssh + 1)
      do jssh = issh, nssh
        kssh = kssh + 1
        if (ilssh /= lsshxc(in1, jssh)) iszero(kssh) = .true.
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
        xnocc_in(issh) = xnocc(issh, in1) + ddq
        xmatt(issh + 1, ix) = ddq
      end do

      ! Obtain the density and respective derivatives needed for evaluating the
      ! exchange-correlation interactions (LDA or GGA).
      ! We have to avoid change xnocc_in !!
      rho1c = 0.0_wp
      rhop1c = 0.0_wp
      rhopp1c = 0.0_wp
      call rho1c_store (in1, nsh_max, nssh, 0.0_wp, 1, drho, rcutoff, &
      &                 xnocc_in, 1, wfmax_points, rho1c, rhop1c, rhopp1c)

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
        call get_potxc1c (iexc, fraction, rho, rh, rhp, rhpp, exc, vxc, &
        &                 dnuxc, dnuxcs, dexcc)
        vxc = hartree*vxc
        exc = hartree*exc
        rho = rho*abohr

        kssh = 0
        do issh = 1, nssh
          do jssh = issh, nssh
            kssh = kssh + 1
            if (iszero(kssh)) cycle
            tmp = psiofr(in1, issh, rho)*psiofr(in1, jssh, rho)*factor*rho*rho
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
    call math_lstsq(xmatt, tmpmat, is_a_trans=.true., info=info)
    deallocate(xmatt)
    if (info /= 0) then
      deallocate(tmpmat, iszero)
      errmsg = 'INTEGRALS/onecenterxc.f90: call to `(x)gelst` failed with info = '
      write (stderr, '(a)') errmsg
      stop info
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
      im = -lsshxc(in1, issh)
      do ilssh = 1, (2*lsshxc(in1, issh) + 1)
        iorb = iorb + 1
        orb2ssh(iorb) = issh
        morbs(iorb) = im
        im = im + 1
      end do
    end do

    ! Write log
    write (auxz,'(i2.2)') nzx(in1)
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
    write (360,'(a)') ''
    do ix = 0, nssh
      write (360,'(100ES15.7)') (exc1crho(ix, issh, issh), issh = 1, nssh)
      write (360,'(a)') ''
    end do
    close(360)
    deallocate(orb2ssh, morbs, exc1crho, nuxc1crho)
  end do ! do in1 = 1, nspec
end subroutine onecenterxc
