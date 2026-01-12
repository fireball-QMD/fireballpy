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
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
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
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! onecenterxc.f90
! Program Description
! ===========================================================================
!       This routine calculates the one-center integrals for the exchange-
! correlation interactions.
! ===========================================================================
! Original code from Juergen Fritsch
!
! Code rewritten by:
! C. Roldan Pinero
! ===========================================================================
subroutine onecenterxc (nspec, nspec_max, nsh_max, wfmax_points,      &
&                       iexc, inux, isnux, fraction, nsshxc, lsshxc,  &
&                       rcutoffa_max, xnocc, dqorb, iderorb, what,    &
&                       signature, drr_rho, nzx)
  use iso_fortran_env, only: error_unit
  use constants, only: abohr, hartree
  use precision, only: long
  implicit none

  ! Argument Declaration and Description
  ! ===========================================================================
  integer, intent (in) :: iexc, nsh_max, nspec, nspec_max, wfmax_points, &
  & inux, isnux
  real(kind=long), intent (in) :: fraction
  character (len=70), intent(in) :: signature
  integer, intent (in), dimension (nspec_max) :: iderorb, nsshxc, nzx
  integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
  real(kind=long), intent (in), dimension (nspec_max) :: dqorb, drr_rho, &
  & rcutoffa_max
  real(kind=long), intent (in), dimension (nsh_max, nspec_max) :: xnocc
  character (len=70), intent (in), dimension (nspec_max) :: what

  ! Local Variable Declaration and Description
  ! ===========================================================================
  integer :: irho, issh, jssh, kssh, in1, lssh, nssh, nssh1, ix, ndq, &
  & nnz, npts, nnrho, nrhs, lwork, info, inuxs
  real(kind=long) :: dnuxc, dnuxcs, drho, exc, dexcc, factor, rcutoff, rho, &
  & rhomax, rhomin, rh, rhp, rhpp, vxc, ddq, tmp
  character(2) :: auxz
  real(kind=long), dimension(nsh_max) :: xnocc_in
  real(kind=long), dimension (:), allocatable :: rho1c, rhop1c, rhopp1c, &
  & work
  real(kind=long), dimension (:,:), allocatable :: xmatt, eexc, vvxc, ddnuxc, &
  & ddnuxcs, tmpmat, fite, fitv, dq
  real(kind=long), dimension(:,:,:), allocatable :: exc1crho, nuxc1crho, &
  & nuxc_onecenter, nuxcs_onecenter
  logical, dimension(:), allocatable :: iszero
  real(kind=long), external :: psiofr
  character(256) :: errmsg

  ! Procedure
  ! ===========================================================================
  ! Allocate known size arrays
  allocate (rho1c (wfmax_points))
  allocate (rhop1c (wfmax_points))
  allocate (rhopp1c (wfmax_points))

  ! There is one additional condition for nuxcs which we need to ensure...
  inuxs = 0
  if (isnux == 1) then
    if (iexc == 11) then
      inuxs = 1
    else
      write (*,*) ''
      write (*,*) 'The spin-polarization option is currently'
      write (*,*) 'implemented only for the LSDA-VWN exchange'
      write (*,*) 'correlation model.'
      write (*,*) 'You have iexc = ', iexc, ' but iexc can only'
      write (*,*) 'equal iexc = 11 for spin-polarization!'
    end if
  end if

  ! Set up the header for the output file.
  open (unit=36, file='coutput/xc1c_dqi.dat', status='unknown')
  open (unit=37, file='coutput/exc1crho.dat', status='unknown')
  open (unit=38, file='coutput/nuxc1crho.dat', status='unknown')
  write (36,100)
  write (36,*) 'All one center matrix elements'
  write (36,*) 'created by:'
  write (36,200) signature
  write (37,100)
  write (37,*) 'All one center matrix elements'
  write (37,*) 'created by:'
  write (37,200) signature
  write (38,100)
  write (38,*) 'All one center matrix elements'
  write (38,*) 'created by:'
  write (38,200) signature
  if (inux == 1) then
    open (unit=39, file='coutput/nuxc_onecenter.dat', status='unknown')
    write (39,100)
    write (39,*) 'All one center matrix elements'
    write (39,*) 'created by:'
    write (39,200) signature
  end if
  if (inuxs == 1) then
    open (unit=40, file='coutput/nuxcs_onecenter.dat', status='unknown')
    write (40,100)
    write (40,*) 'All one center matrix elements (spin)'
    write (40,*) 'created by:'
    write (40,200) signature
  end if
  do in1 = 1, nspec
    write (36,300) what(in1)
    write (37,300) what(in1)
    write (38,300) what(in1)
    write (auxz,'(i2.2)') nzx(in1)
    open (unit=360, file='coutput/xc1c_dqi.'//auxz//'.dat', status='unknown')
    open (unit=370, file='coutput/exc1crho.'//auxz//'.dat', status='unknown')
    open (unit=380, file='coutput/nuxc1crho.'//auxz//'.dat', status='unknown')
    write (360,100)
    write (360,*) 'Z = ', nzx(in1), ' one center matrix elements'
    write (360,*) 'created by:'
    write (360,200) signature
    write (360,300) what(in1)
    write (360,100)
    close (360)
    write (370,100)
    write (370,*) 'Z = ', nzx(in1), ' one center matrix elements'
    write (370,*) 'created by:'
    write (370,200) signature
    write (370,300) what(in1)
    write (370,100)
    close (370)
    write (380,100)
    write (380,*) 'Z = ', nzx(in1), ' one center matrix elements'
    write (380,*) 'created by:'
    write (380,200) signature
    write (380,300) what(in1)
    write (380,100)
    close (380)
    if (inux == 1) then
      open (unit=390, file='coutput/nuxc_onecenter.'//auxz//'.dat', status='unknown')
      write (390,*) 'Z = ', nzx(in1), ' one center matrix elements'
      write (390,*) 'created by:'
      write (390,200) signature
      write (390,300) what(in1)
      write (390,100)
      close (390)
    end if
    if (inuxs == 1) then
      open (unit=400, file='coutput/nuxcs_onecenter.'//auxz//'.dat', status='unknown')
      write (400,*) 'Z = ', nzx(in1), ' one center matrix elements (spin)'
      write (400,*) 'created by:'
      write (400,200) signature
      write (400,300) what(in1)
      write (400,100)
      close (400)
    end if
  end do
  write (36,100)
  write (37,100)
  write (38,100)
  close (36)
  close (37)
  close (38)
  if (inux == 1) then
    write (39,100)
    close (39)
  end if
  if (inuxs == 1) then
    write (40,100)
    close (40)
  end if

  do in1 = 1, nspec
    nssh = nsshxc(in1)
    nssh1 = nssh + 1

    ! Needed for charge corrections
    drho = drr_rho(in1)
    rcutoff = rcutoffa_max(in1)
    rhomin = 0.0_long
    rhomax = rcutoff
    nnrho = nint((rhomax - rhomin)/drho) + 1

    ! Prepare the increments
    ndq = 3
    npts = ndq**nssh
    nnz = (nssh*(nssh + 1))/2
    allocate(dq(nssh, ndq))
    do issh = 1,nssh
      lssh = lsshxc(in1,issh)
      if (lssh == 0) then
        if (issh == 1) then
          dq(issh,1) = -0.10_long
          dq(issh,2) = 0.00_long
          dq(issh,3) = 0.10_long
        else
          dq(issh,1) = 0.00_long
          dq(issh,2) = 0.05_long
          dq(issh,3) = 0.10_long
        end if
      else if (lssh == 1) then
        if (issh == 2) then
          dq(issh,1) = -0.20_long
          dq(issh,2) = 0.00_long
          dq(issh,3) = 0.20_long
        else
          dq(issh,1) = 0.00_long
          dq(issh,2) = 0.10_long
          dq(issh,3) = 0.20_long
        end if
      else
        if (issh == 3) then
          dq(issh,1) = -0.50_long
          dq(issh,2) = 0.00_long
          dq(issh,3) = 0.50_long
        else
          dq(issh,1) = 0.00_long
          dq(issh,2) = 0.25_long
          dq(issh,3) = 0.50_long
        end if
      end if
    end do

    ! We get the full grid of values first
    allocate(iszero(nnz))
    iszero = .false.
    kssh = 0
    do issh = 1, nssh
      do jssh = issh, nssh
        kssh = kssh + 1
        if (lsshxc(in1, issh) /= lsshxc(in1, jssh)) iszero(kssh) = .true.
      end do
    end do
    allocate(eexc(npts, nnz), vvxc(npts, nnz))
    eexc = 0.0_long
    vvxc = 0.0_long
    if (inux == 1) then
      allocate(ddnuxc(npts, nnz))
      ddnuxc = 0.0_long
    end if
    if (inuxs == 1) then
      allocate(ddnuxcs(npts, nnz))
      ddnuxcs = 0.0_long
    end if
    allocate(xmatt(nssh1, npts))
    do ix = 1, npts
      ! Set charges
      xmatt(1, ix) = 1.0_long
      do issh = 1, nssh
        ddq = dq(issh, 1 + mod(ix - 1, ndq**issh)/ndq**(issh - 1))
        xnocc_in(issh) = xnocc(issh, in1) + ddq
        xmatt(issh + 1, ix) = ddq
      end do

      ! Obtain the density and respective derivatives needed for evaluating the
      ! exchange-correlation interactions (LDA or GGA).
      ! We have to avoid change xnocc_in !!
      rho1c = 0.0_long
      rhop1c = 0.0_long
      rhopp1c = 0.0_long
      call rho1c_store (in1, nsh_max, nssh, 0.0_long, 1, drho, rcutoff, &
      &                 xnocc_in, 1, wfmax_points, rho1c, rhop1c, rhopp1c)

      ! Integrals <i|exc(i)|i> and <i.nu|mu(i)|i.nu'>
      do irho = 1, nnrho
        rho = rhomin + real(irho - 1, kind=long)*drho
        factor = 0.66666666666666666667_long*drho
        if (mod(irho, 2) == 0) factor = 2.0_long*factor
        if (irho == 1 .or. irho == nnrho) factor = 0.5_long*factor

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

        ! Important note. We have ONLY programmed in ldaxc.f the nu potential
        ! (nu = dmu/dn) for iexc = 3. So even if you do a different LDA, or even GGA,
        ! we will use Ceperly/Alder (iexc = 3) for the nu part.
        if (inux == 1) then
          if (iexc /= 1 .and. iexc /= 11) then
            rho = rho/abohr
            call get_potxc1c (3, fraction, rho, rh, rhp, rhpp, exc, vxc, &
            & dnuxc, dnuxcs, dexcc)
            rho = rho*abohr
          end if
          dnuxc = hartree*dnuxc*abohr**3

          kssh = 0
          do issh = 1, nssh
            do jssh = issh, nssh
              kssh = kssh + 1
              if (iszero(kssh)) cycle
              tmp = psiofr(in1, issh, rho)*psiofr(in1, jssh, rho)*factor*rho*rho
              ddnuxc(ix, kssh) = ddnuxc(ix, kssh) + tmp*dnuxc
            end do
          end do
        end if

        if (inuxs == 1) then
          dnuxcs = hartree*dnuxcs*abohr**3
          kssh = 0
          do issh = 1, nssh
            do jssh = issh, nssh
              kssh = kssh + 1
              if (iszero(kssh)) cycle
              tmp = psiofr(in1, issh, rho)*psiofr(in1, jssh, rho)*factor*rho*rho
              ddnuxcs(ix, kssh) = ddnuxcs(ix, kssh) + tmp*dnuxcs
            end do
          end do
        end if
      end do
    end do
    deallocate(dq)

    ! Perform the fits
    nrhs = 2
    if (inux == 1) nrhs = nrhs + 1
    if (inuxs == 1) nrhs = nrhs + 1
    allocate(tmpmat(npts, nrhs*nnz))
    tmpmat(:, 1:nnz) = eexc
    tmpmat(:, nnz+1:2*nnz) = vvxc
    deallocate(eexc, vvxc)
    nrhs = 2
    if (inux == 1) then
      tmpmat(:, nrhs*nnz+1:(nrhs+1)*nnz) = ddnuxc
      nrhs = nrhs + 1
      deallocate(ddnuxc)
    end if
    if (inuxs == 1) then
      tmpmat(:, nrhs*nnz+1:(nrhs+1)*nnz) = ddnuxcs
      nrhs = nrhs + 1
      deallocate(ddnuxcs)
    end if

    ! Query does not work, don't ask why
    lwork = nssh1 + max(nssh1, nrhs*nnz)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate(work(lwork))
    call dgelst('T', nssh1, npts, nrhs*nnz, xmatt, nssh1, tmpmat, npts, work, lwork, info)
    deallocate(work, xmatt)
    if (info /= 0) then
      deallocate(tmpmat, iszero, rho1c, rhop1c, rhopp1c)
      errmsg = 'INTEGRALS/onecenterxc.f90: call to `dgelst` failed with info = '
      write (error_unit, '(a)') errmsg
      stop info
    end if

    ! Prepare for output
    allocate(exc1crho(0:nssh,nssh,nssh), nuxc1crho(0:nssh,nssh,nssh))
    exc1crho = 0.0_long
    nuxc1crho = 0.0_long
    do ix = 0, nssh
      kssh = 0
      do issh = 1, nssh
        do jssh = issh, nssh
          kssh = kssh + 1
          if (iszero(kssh)) cycle
          exc1crho(ix, issh, jssh) = tmpmat(ix+1, kssh)
          exc1crho(ix, jssh, issh) = tmpmat(ix+1, kssh)
          nuxc1crho(ix, issh, jssh) = tmpmat(ix+1, kssh+nnz)
          nuxc1crho(ix, jssh, issh) = tmpmat(ix+1, kssh+nnz)
        end do
      end do
    end do
    nrhs = 2
    if (inux == 1) then
      allocate(nuxc_onecenter(0:nssh,nssh,nssh))
      do ix = 0, nssh
        kssh = 0
        do issh = 1, nssh
          do jssh = issh, nssh
            kssh = kssh + 1
            if (iszero(kssh)) cycle
            nuxc_onecenter(ix, issh, jssh) = tmpmat(ix+1, kssh+nrhs*nnz)
            nuxc_onecenter(ix, jssh, issh) = tmpmat(ix+1, kssh+nrhs*nnz)
          end do
        end do
      end do
      nrhs = nrhs + 1
    end if
    if (inuxs == 1) then
      allocate(nuxcs_onecenter(0:nssh,nssh,nssh))
      do ix = 0, nssh
        kssh = 0
        do issh = 1, nssh
          do jssh = issh, nssh
            kssh = kssh + 1
            if (iszero(kssh)) cycle
            nuxcs_onecenter(ix, issh, jssh) = tmpmat(ix+1, kssh+nrhs*nnz)
            nuxcs_onecenter(ix, jssh, issh) = tmpmat(ix+1, kssh+nrhs*nnz)
          end do
        end do
      end do
    end if
    deallocate(tmpmat, iszero)

    ! Write log
    write (auxz,'(i2.2)') nzx(in1)
    write (*,*) ' '
    write (*,*) 'Writing output to: coutput/xc1c_dqi.'//auxz//'.dat'
    write (*,*) 'Writing output to: coutput/exc1crho.'//auxz//'.dat'
    write (*,*) 'Writing output to: coutput/nuxc1crho.'//auxz//'.dat'
    if (inux == 1) write (*,*) 'Writing output to: coutput/nuxc_onecenter.'//auxz//'.dat'
    if (inuxs == 1) write (*,*) 'Writing output to: coutput/nuxcs_onecenter.'//auxz//'.dat'
    write (*,*) ' '

    ! Write output
    open (unit=36, file='coutput/xc1c_dqi.dat', position='append', status='old')
    open (unit=360, file='coutput/xc1c_dqi.'//auxz//'.dat', position='append', status='old')
    write (36,400) in1, nssh
    write (360,400) in1, nssh
    do issh = 1, nssh
      write (36,501) (exc1crho(0,issh,jssh), jssh = 1, nssh)
      write (360,501) (exc1crho(0,issh,jssh), jssh = 1, nssh)
    end do
    write (36,*)
    write (360,*)
    do issh = 1, nssh
      write (36,501) (nuxc1crho(0,issh,jssh), jssh = 1, nssh)
      write (360,501) (nuxc1crho(0,issh,jssh), jssh = 1, nssh)
    end do
    close(36)
    close(360)

    open (unit=37, file='coutput/exc1crho.dat', position='append', status='old')
    open (unit=38, file='coutput/nuxc1crho.dat', position='append', status='old')
    open (unit=370, file='coutput/exc1crho.'//auxz//'.dat', position='append', status='old')
    open (unit=380, file='coutput/nuxc1crho.'//auxz//'.dat', position='append', status='old')
    do kssh = 1, nssh
      write (37,410) in1, nssh, kssh
      write (370,410) in1, nssh, kssh
      write (38,410) in1, nssh, kssh
      write (380,410) in1, nssh, kssh
      do issh = 1, nssh
        write (37,501) (exc1crho(kssh,issh,jssh), jssh = 1, nssh)
        write (370,501) (exc1crho(kssh,issh,jssh), jssh = 1, nssh)
        write (38,501) (nuxc1crho(kssh,issh,jssh), jssh = 1, nssh)
        write (380,501) (nuxc1crho(kssh,issh,jssh), jssh = 1, nssh)
      end do
    end do
    close(37)
    close(38)
    close(370)
    close(380)
    deallocate(exc1crho, nuxc1crho)

    if (inux == 1) then
      open (unit=39, file='coutput/nuxc_onecenter.dat', position='append', status='old')
      open (unit=390, file='coutput/nuxc_onecenter.'//auxz//'.dat', position='append', status='old')
      do kssh = 1, nssh
        write (39,410) in1, nssh, kssh
        write (390,410) in1, nssh, kssh
        do issh = 1, nssh
          write (39,501) (nuxc_onecenter(kssh,issh,jssh), jssh = 1, nssh)
          write (390,501) (nuxc_onecenter(kssh,issh,jssh), jssh = 1, nssh)
        end do
      end do
      close(39)
      close(390)
      deallocate(nuxc_onecenter)
    end if
    if (inuxs == 1) then
      open (unit=40, file='coutput/nuxcs_onecenter.dat', position='append', status='old')
      open (unit=400, file='coutput/nuxcs_onecenter.'//auxz//'.dat', &
      &  position='append', status='old')
      do kssh = 1, nssh
        write (40,410) in1, nssh, kssh
        write (400,410) in1, nssh, kssh
        do issh = 1, nssh
          write (40,501) (nuxcs_onecenter(kssh,issh,jssh), jssh = 1, nssh)
          write (400,501) (nuxcs_onecenter(kssh,issh,jssh), jssh = 1, nssh)
        end do
      end do
      close(40)
      close(400)
      deallocate(nuxcs_onecenter)
    end if
  end do ! do in1 = 1, nspec
  deallocate (rho1c, rhop1c, rhopp1c)

  ! Format Statements
  ! ===========================================================================
100 format (70('='))
200 format (2x, a45)
300 format (a70)
400 format (2x, i3, 2x, i3)
410 format (2x, i3, 2x, i3, 2x, i3)
450 format (2x, i3, 2x, i3, 2x, i3)
470 format (2x, f7.3, 1x, f7.3, 1x, f7.3)
480 format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
500 format (8e20.10)
501 format (8e20.10)
600 format (1x, i3, 2x, f10.5)

end subroutine onecenterxc
