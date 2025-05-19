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
  &                     iexc, fraction, nsshxc, lsshxc, rcutoffa_max, &
  &                     xnocc, dqorb, iderorb, what, signature,       &
  &                     drr_rho, nzx)
  use constants
  use precision
  implicit none

  ! Argument Declaration and Description
  ! ===========================================================================
  integer, intent (in) :: iexc, nsh_max, nspec, nspec_max, wfmax_points
  real(kind=long), intent (in) :: fraction
  character (len=70) :: signature
  integer, intent (in), dimension (nspec_max) :: iderorb, nsshxc, nzx
  integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
  real(kind=long), intent (in), dimension (nspec_max) :: dqorb, drr_rho, &
    & rcutoffa_max
  real(kind=long), intent (in), dimension (nsh_max, nspec_max) :: xnocc
  character (len=70), intent (in), dimension (nspec_max) :: what

  ! Local Variable Declaration and Description
  ! ===========================================================================
  integer :: irho, issh, jssh, kssh, in1, lssh, nssh, nssh1, ix, ndq, &
    & nnz, npts, nnrho, lwork, info
  real(kind=long) :: dnuxc, dnuxcs, drho, exc, dexcc, factor, rcutoff, rho, &
    & rhomax, rhomin, rh, rhp, rhpp, vxc, ddq, tmp
  character(2) :: auxz
  integer, dimension(:), allocatable :: ipiv
  real(kind=long), dimension(nsh_max) :: xnocc_in
  real(kind=long), dimension (:), allocatable :: rho1c, rhop1c, rhopp1c, &
    & work
  real(kind=long), dimension (:,:), allocatable :: xmatt, eexc, vvxc, tmpmat, &
    & fite, fitv, dq
  real(kind=long), dimension(:,:,:), allocatable :: exc1crho, nuxc1crho
  logical, dimension(:), allocatable :: iszero
  real(kind=long), external :: psiofr

  ! Procedure
  ! ===========================================================================
  ! Allocate known size arrays
  allocate (rho1c (wfmax_points))
  allocate (rhop1c (wfmax_points))
  allocate (rhopp1c (wfmax_points))

  ! Set up the header for the output file.
  open (unit=36, file='coutput/xc1c_dqi.dat', status='unknown')
  open (unit=37, file='coutput/exc1crho.dat', status='unknown')
  open (unit=38, file='coutput/nuxc1crho.dat', status='unknown')
  write (36,100)
  write (36,*) ' All one center matrix elements '
  write (36,*) ' created by: '
  write (36,200) signature
  write (37,100)
  write (37,*) ' All one center matrix elements '
  write (37,*) ' created by: '
  write (37,200) signature
  write (38,100)
  write (38,*) ' All one center matrix elements '
  write (38,*) ' created by: '
  write (38,200) signature
  do in1 = 1, nspec
    write (36,300) what(in1)
    write (37,300) what(in1)
    write (38,300) what(in1)
    write (auxz,'(i2.2)') nzx(in1)
    open (unit=360, file='coutput/xc1c_dqi.'//auxz//'.dat', status='unknown')
    open (unit=370, file='coutput/exc1crho.'//auxz//'.dat', status='unknown')
    open (unit=380, file='coutput/nuxc1crho.'//auxz//'.dat', status='unknown')
    write (360,100)
    write (360,*) ' Z = ', nzx(in1), ' one center matrix elements '
    write (360,*) ' created by: '
    write (360,200) signature
    write (360,300) what(in1)
    write (360,100)
    close (360)
    write (370,100)
    write (370,*) ' Z = ', nzx(in1), ' one center matrix elements '
    write (370,*) ' created by: '
    write (370,200) signature
    write (370,300) what(in1)
    write (370,100)
    close (370)
    write (380,100)
    write (380,*) ' Z = ', nzx(in1), ' one center matrix elements '
    write (380,*) ' created by: '
    write (380,200) signature
    write (380,300) what(in1)
    write (380,100)
    close (380)
  end do
  write (36,100)
  write (37,100)
  write (38,100)
  close (36)
  close (37)
  close (38)

  do in1 = 1, nspec
    nssh = nsshxc(in1)
    nssh1 = nssh + 1

    ! Needed for charge corrections
    drho = drr_rho(in1)
    rcutoff = rcutoffa_max(in1)
    rhomin = 0.0d0
    rhomax = rcutoff
    nnrho = nint((rhomax - rhomin)/drho) + 1

    ! Prepare the increments
    ndq = 3
    npts = ndq**nssh
    nnz = (nssh*(nssh + 1))/2
    allocate(dq(nssh, ndq))
    do issh = 1,nssh
      lssh = lsshxc(in1,issh)
      if (lssh .eq. 0) then
        if (issh .eq. 1) then
          dq(issh,1) = -0.10d0
          dq(issh,2) = 0.00d0
          dq(issh,3) = 0.10d0
        else
          dq(issh,1) = 0.00d0
          dq(issh,2) = 0.05d0
          dq(issh,3) = 0.10d0
        end if
      else if (lssh .eq. 1) then
        if (issh .eq. 2) then
          dq(issh,1) = -0.20d0
          dq(issh,2) = 0.00d0
          dq(issh,3) = 0.20d0
        else
          dq(issh,1) = 0.00d0
          dq(issh,2) = 0.10d0
          dq(issh,3) = 0.20d0
        end if
      else
        if (issh .eq. 3) then
          dq(issh,1) = -0.50d0
          dq(issh,2) = 0.00d0
          dq(issh,3) = 0.50d0
        else
          dq(issh,1) = 0.00d0
          dq(issh,2) = 0.25d0
          dq(issh,3) = 0.50d0
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
        if (lsshxc(in1, issh) .ne. lsshxc(in1, jssh)) iszero(kssh) = .true.
      end do
    end do
    allocate(eexc(npts, nnz), vvxc(npts, nnz))
    allocate(xmatt(nssh1, npts))
    eexc = 0.0d0
    vvxc = 0.0d0
    do ix = 1, npts
      ! Set charges
      xmatt(1, ix) = 1.0d0
      do issh = 1, nssh
        ddq = dq(issh, 1 + mod(ix - 1, ndq**issh)/ndq**(issh - 1))
        xnocc_in(issh) = xnocc(issh, in1) + ddq
        xmatt(issh + 1, ix) = ddq
      end do

      ! Obtain the density and respective derivatives needed for evaluating the
      ! exchange-correlation interactions (LDA or GGA).
      ! We have to avoid change xnocc_in !!
      rho1c = 0.0d0
      rhop1c = 0.0d0
      rhopp1c = 0.0d0
      call rho1c_store (in1, nsh_max, nssh, 0.0d0, 1, drho, rcutoff, &
        & xnocc_in, 1, wfmax_points, rho1c, rhop1c, rhopp1c)

      ! Integrals <i|exc(i)|i> and <i.nu|mu(i)|i.nu'>
      do irho = 1, nnrho
        rho = rhomin + real(irho - 1, kind=long)*drho
        factor = 0.66666666666666666667d0*drho
        if (mod(irho, 2) .eq. 0) factor = 2.0d0*factor
        if (irho .eq. 1 .or. irho .eq. nnrho) factor = 0.5d0*factor

        ! Compute the exchange correlation potential
        rho = rho/abohr
        rh = rho1c(irho)*abohr**3
        rhp = rhop1c(irho)*abohr**4
        rhpp = rhopp1c(irho)*abohr**5
        call get_potxc1c (iexc, fraction, rho, rh, rhp, rhpp, exc, vxc, &
          & dnuxc, dnuxcs, dexcc)
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

    ! Perform the fit. The (XX**T)X may be reutilised
    allocate(ipiv(nssh1))
    allocate(tmpmat(nssh1, nssh1))
    call dgemm('N', 'T', nssh1, nssh1, npts, 1.0d0, xmatt, nssh1, &
      & xmatt, nssh1, 0.0d0, tmpmat, nssh1)
    call dgetrf(nssh1, nssh1, tmpmat, nssh1, ipiv, info)
    if (info .ne. 0) then
      write (*,*) 'ERROR in onecenterxc.f90'
    end if
    allocate(work(1))
    call dgetri(nssh1, tmpmat, nssh1, ipiv, work, -1, info)
    lwork = nint(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dgetri(nssh1, tmpmat, nssh1, ipiv, work, lwork, info)
    if (info .ne. 0) then
      write (*,*) 'ERROR in onecenterxc.f90'
    end if
    call dgemm('N', 'N', nssh1, npts, nssh1, 1.0d0, tmpmat, nssh1, &
      & xmatt, nssh1, 0.0d0, xmatt, nssh1)
    deallocate(ipiv, work, tmpmat)

    allocate(fite(nssh1, nnz), fitv(nssh1, nnz))
    do kssh = 1, nnz
      if (iszero(kssh)) then
        fite(:, kssh) = 0.0d0
        fitv(:, kssh) = 0.0d0
        cycle
      end if
      call dgemv('N', nssh1, npts, 1.0d0, xmatt, nssh1, &
        & eexc(:, kssh), 1, 0.0d0, fite(:, kssh), 1)
      call dgemv('N', nssh1, npts, 1.0d0, xmatt, nssh1, &
        & vvxc(:, kssh), 1, 0.0d0, fitv(:, kssh), 1)
    end do
    deallocate(eexc, vvxc, xmatt)

    ! Prepare for output
    allocate(exc1crho(nssh1,nssh,nssh), nuxc1crho(nssh1,nssh,nssh))
    exc1crho = 0.0d0
    nuxc1crho = 0.0d0
    do ix = 1, nssh1
      kssh = 0
      do issh = 1, nssh
        do jssh = issh, nssh
          kssh = kssh + 1
          if (iszero(kssh)) cycle
          exc1crho(ix, issh, jssh) = fite(ix, kssh)
          exc1crho(ix, jssh, issh) = fite(ix, kssh)
          nuxc1crho(ix, issh, jssh) = fitv(ix, kssh)
          nuxc1crho(ix, jssh, issh) = fitv(ix, kssh)
        end do
      end do
    end do
    deallocate(iszero, fite, fitv)

    ! Write log
    write (auxz,'(i2.2)') nzx(in1)
    write (*,*) '  '
    write (*,*) ' Writing output to: coutput/xc1c_dqi.'//auxz//'.dat '
    write (*,*) ' Writing output to: coutput/exc1crho.'//auxz//'.dat '
    write (*,*) ' Writing output to: coutput/nuxc1crho.'//auxz//'.dat '
    write (*,*) '  '

    ! Write output
    open (unit=36, file='coutput/xc1c_dqi.dat', position='append', status='old')
    open (unit=360, file='coutput/xc1c_dqi.'//auxz//'.dat', position='append', status='old')
    write (36,400) in1, nssh
    write (360,400) in1, nssh
    do issh = 1, nssh
      write (36,501) (exc1crho(1,issh,jssh), jssh = 1, nssh)
      write (360,501) (exc1crho(1,issh,jssh), jssh = 1, nssh)
    end do
    write (36,*)
    write (360,*)
    do issh = 1, nssh
      write (36,501) (nuxc1crho(1,issh,jssh), jssh = 1, nssh)
      write (360,501) (nuxc1crho(1,issh,jssh), jssh = 1, nssh)
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
        write (37,501) (exc1crho(kssh+1,issh,jssh), jssh = 1, nssh)
        write (370,501) (exc1crho(kssh+1,issh,jssh), jssh = 1, nssh)
        write (38,501) (nuxc1crho(kssh+1,issh,jssh), jssh = 1, nssh)
        write (380,501) (nuxc1crho(kssh+1,issh,jssh), jssh = 1, nssh)
      end do
    end do
    close(37)
    close(38)
    close(370)
    close(380)
    deallocate(exc1crho, nuxc1crho)
  end do ! do in1 = 1, nspec
  deallocate (rho1c, rhop1c, rhopp1c)

  ! Format Statements
  ! ===========================================================================
  100     format (70('='))
  200     format (2x, a45)
  300     format (a70)
  400     format (2x, i3, 2x, i3)
  410     format (2x, i3, 2x, i3, 2x, i3)
  450     format (2x, i3, 2x, i3, 2x, i3)
  470     format (2x, f7.3, 1x, f7.3, 1x, f7.3)
  480     format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
  500     format (8d20.10)
  501     format (8d20.10)
  600     format (1x, i3, 2x, f10.5)

  return
end subroutine onecenterxc
