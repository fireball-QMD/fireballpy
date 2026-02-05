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
! University of Utah - James P. Lewis Kurt R. Glaesemann
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


! twocenterxc.f90
! Program Description
! ===========================================================================
!      This code is a general two-center integration routine for matrix
! elements of the form <psi(1)|V(1)|psi(2)>.  Thus, V(1) is located at the
! site of one of the orbitals.  The potential V(1) is something like Vxc for
! the exchange correlation potential, Vna for the neutral atom potential, or
! 1 for the overlap term.
!
! ===========================================================================
! Code written by C. Roldán Piñero
! ===========================================================================
module twocenterxc
  use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
  use precision, only: wp
  use utils, only: utils_progress_bar, utils_clean_progress_bar
  implicit none
  private

  public :: TWOCENTER_XC
  public :: twocenter_init
  public :: twocenter_calc

  integer, parameter :: TWOCENTER_XC = 1    ! 2^0
  integer, parameter :: TWOCENTER_NUM_INTERACTION = 1
  logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTION)

  integer :: iexc_, ispec1_, ispec2_, nzx1_, nzx2_, nz_, nrho_, nnrho_, ndd_,  &
  &          index_max_
  real(kind=wp) :: fraction_, rcutoff1_, rcutoff2_, rhomin_, rhomax_, dz_, drho_
  character(len=2) :: atom1_, atom2_
  character(len=70) :: what1_, what2_, signature_
  integer, allocatable :: nleft_(:), lleft_(:), mleft_(:), nright_(:),         &
  &                       lright_(:), mright_(:)

contains

  ! TODO: delete non-sense and they should come from a module
  integer function twocenter_init(iexc, fraction, ispec1, ispec2, atom1,       &
  &                               atom2, what1, what2, nzx1, nzx2, rcutoff1,   &
  &                               rcutoff2, nz, nrho, ndd,                     &
  &                               index_max, nleft, lleft, mleft, nright,      &
  &                               lright, mright, signature)
    use precision, only: wp
    integer, intent(in) :: iexc, ispec1, ispec2, nzx1, nzx2, nz,               &
    &                      nrho, ndd, index_max, nleft(:), lleft(:), mleft(:), &
    &                      nright(:), lright(:), mright(:)
    real(kind=wp), intent(in) :: fraction, rcutoff1, rcutoff2
    character(len=2), intent(in) :: atom1, atom2
    character(len=70), intent(in) :: what1, what2, signature

    if (ndd < 107) then
      write (stderr, '(a)') '[ERROR] twocenterxc.f90: ndd must be at least 107'
      twocenter_init = -15
      return
    end if

    if (allocated(nleft_)) deallocate(nleft_)
    if (allocated(lleft_)) deallocate(lleft_)
    if (allocated(mleft_)) deallocate(mleft_)
    if (allocated(nright_)) deallocate(nright_)
    if (allocated(lright_)) deallocate(lright_)
    if (allocated(mright_)) deallocate(mright_)

    allocate(nleft_(index_max), lleft_(index_max), mleft_(index_max),          &
    & nright_(index_max), lright_(index_max), mright_(index_max),              &
    & stat=twocenter_init)
    if (twocenter_init /= 0) return

    iexc_ = iexc
    ispec1_ = ispec1
    ispec2_ = ispec2
    nzx1_ = nzx1
    nzx2_ = nzx2
    nz_ = nz
    nrho_ = nrho
    ndd_ = ndd
    index_max_ = index_max
    fraction_ = fraction
    rcutoff1_ = rcutoff1
    rcutoff2_ = rcutoff2
    atom1_ = atom1
    atom2_ = atom2
    what1_ = what1
    what2_ = what2
    signature_ = signature
    nleft_ = nleft(1:index_max_)
    lleft_ = lleft(1:index_max_)
    mleft_ = mleft(1:index_max_)
    nright_ = nright(1:index_max_)
    lright_ = lright(1:index_max_)
    mright_ = mright(1:index_max_)
    rhomin_ = 0.0_wp
    rhomax_ = min(rcutoff1_, rcutoff2_)
    dz_ = 0.5_wp*(rcutoff1_ + rcutoff2_)/real(nz, kind=wp)
    drho_ = max(rcutoff1_, rcutoff2_)/real(nrho, kind=wp)
    nnrho_ = int((rhomax_ - rhomin_)/drho_)
    if (iand(nnrho_, 1) == 0) nnrho_ = nnrho_ + 1

    twocenter_init = 0
    return
  end function twocenter_init

  subroutine twocenter_get_interactions(interactions)
    integer, intent(in) :: interactions
    integer :: i, p, inter
    twocenter_interactions = .false.
    inter = interactions
    p = ishft(1, TWOCENTER_NUM_INTERACTION - 1)
    do i = TWOCENTER_NUM_INTERACTION, 1, -1
      if (mod(inter, p) == 0) then
        twocenter_interactions(i) = .true.
        inter = inter - p
      end if
      p = ishft(p, -1)
    end do
  end subroutine twocenter_get_interactions

  ! TODO: this should come from somewhere else
  subroutine twocenter_charge_grid(dq1, dq2)
    real(kind=wp), intent(out) :: dq1(:,:), dq2(:,:)
    include '../parameters.inc'
    include '../exchange.inc'
    logical :: first_orb(0:3)
    integer :: issh, ilssh, ndq
    ndq = 2

    ! ispec1
    first_orb = .true.
    do issh = 1, nsshxc(ispec1_)
      ilssh = lsshxc(ispec1_, issh)
      if (ilssh == 0) then
        if (first_orb(ilssh)) then
          dq1(issh, 1) = -0.10_wp
          ! dq1(issh, 2) = 0.00_wp
          dq1(issh, 2) = 0.10_wp
          first_orb(ilssh) = .false.
        else
          dq1(issh, 1) = 0.00_wp
          ! dq1(issh, 2) = 0.05_wp
          dq1(issh, 2) = 0.10_wp
        end if
      else if (ilssh == 1) then
        if (first_orb(ilssh)) then
          dq1(issh, 1) = -0.20_wp
          ! dq1(issh, 2) = 0.00_wp
          dq1(issh, 2) = 0.20_wp
          first_orb(ilssh) = .false.
        else
          dq1(issh, 1) = 0.00_wp
          ! dq1(issh, 2) = 0.10_wp
          dq1(issh, 2) = 0.20_wp
        end if
      else
        if (first_orb(ilssh)) then
          dq1(issh, 1) = -0.50_wp
          ! dq1(issh, 2) = 0.00_wp
          dq1(issh, 2) = 0.50_wp
          first_orb(ilssh) = .false.
        else
          dq1(issh, 1) = 0.00_wp
          ! dq1(issh, 2) = 0.25_wp
          dq1(issh, 2) = 0.50_wp
        end if
      end if
    end do

    ! ispec2
    first_orb = .true.
    do issh = 1, nsshxc(ispec2_)
      ilssh = lsshxc(ispec2_, issh)
      if (ilssh == 0) then
        if (first_orb(ilssh)) then
          dq2(issh, 1) = -0.10_wp
          ! dq2(issh, 2) = 0.00_wp
          dq2(issh, 2) = 0.10_wp
          first_orb(ilssh) = .false.
        else
          dq2(issh, 1) = 0.00_wp
          ! dq2(issh, 2) = 0.05_wp
          dq2(issh, 2) = 0.10_wp
        end if
      else if (ilssh == 1) then
        if (first_orb(ilssh)) then
          dq2(issh, 1) = -0.20_wp
          ! dq2(issh, 2) = 0.00_wp
          dq2(issh, 2) = 0.20_wp
          first_orb(ilssh) = .false.
        else
          dq2(issh, 1) = 0.00_wp
          ! dq2(issh, 2) = 0.10_wp
          dq2(issh, 2) = 0.20_wp
        end if
      else
        if (first_orb(ilssh)) then
          dq2(issh, 1) = -0.50_wp
          ! dq2(issh, 2) = 0.00_wp
          dq2(issh, 2) = 0.50_wp
          first_orb(ilssh) = .false.
        else
          dq2(issh, 1) = 0.00_wp
          ! dq2(issh, 2) = 0.25_wp
          dq2(issh, 2) = 0.50_wp
        end if
      end if
    end do
  end subroutine twocenter_charge_grid

  subroutine twocenter_write_header(iounit, dmax, rev)
    integer, intent(in) :: iounit
    real(kind=wp), intent(in) :: dmax
    logical, intent(in) :: rev
    write (iounit, '(a)') ' created by: '
    write (iounit, '(2x,a)') signature_
    if (rev) then
      write (iounit, '(a)') what2_
      write (iounit, '(a)') what1_
    else
      write (iounit, '(a)') what1_
      write (iounit, '(a)') what2_
    end if
    write (iounit, '(a)') repeat('=', 70)
    write (iounit, '(a,f8.4,a,i4,a,i4,a,i4)') '  R(d) = ', dmax,               &
    &     '(A)   ndd = ', ndd_, '   nz_points = ', nz_, '   nrhopoints = ', nrho_
    write (iounit, '(a)') repeat('=', 70)
    if (rev) then
      write (iounit, '(i5,f9.4,a,a)') nzx2_, rcutoff2_, ' <=== ', atom2_
      write (iounit, '(i5,f9.4,a,a)') nzx1_, rcutoff1_, ' <=== ', atom1_
    else
      write (iounit, '(i5,f9.4,a,a)') nzx1_, rcutoff1_, ' <=== ', atom1_
      write (iounit, '(i5,f9.4,a,a)') nzx2_, rcutoff2_, ' <=== ', atom2_
    end if
    write (iounit, '(2x,f8.4,2x,i4)') dmax, ndd_
  end subroutine twocenter_write_header

  integer function twocenter_calc(interactions)
    use math, only: math_lstsq
    integer, intent(in) :: interactions
    include '../parameters.inc'
    include '../exchange.inc'
    include '../pseudopotentials.inc'
    integer :: igrid, index, index_coulomb, issh, l1, l2, m1, m2, n1, n2,      &
    &          io, ix, ix1, ix2, ndq, npts, npts1, npts2, nssh1, nssh2
    real(kind=wp) :: d, dmax, dr, ddq, xnocc_in(nsh_max, 2)
    character(len=40) :: fname1, fname2
    real(kind=wp), allocatable :: dq1(:,:), dq2(:,:), xmatt(:,:), vvxc(:,:)
    character(len=40), allocatable :: fnamel1(:), fnamer1(:), fnamel2(:),      &
    &                                 fnamer2(:)

    call twocenter_get_interactions(interactions)
    if (.not. any(twocenter_interactions)) then
      twocenter_calc = 0
      return
    end if

    nssh1 = nsshxc(ispec1_)
    nssh2 = nsshxc(ispec2_)
    dmax = rcutoff1_ + rcutoff2_
    dr = dmax/real(ndd_ - 1, kind=wp)

    if (twocenter_interactions(1)) then
      allocate(fnamel1(nssh1), fnamel2(nssh2), fnamer1(nssh2), fnamer2(nssh1), &
      &        stat=twocenter_calc)
      if (twocenter_calc /= 0) return
      io = 3610
      write (fname1, '(a,i2.2,a,i2.2,a)') 'vxc_2c.', nzx1_, '.',               &
      &     nzx2_, '.dat'
      write (stdout, '(a)') 'Computing '//trim(fname1)//'...'
      open (unit=io, file='coutput/'//trim(fname1), status='new',              &
      &     action='write', iostat=twocenter_calc)
      if (twocenter_calc /= 0) then
        write (stderr, '(a)') '[ERROR] failed to open '//trim(fname1)
        return
      end if
      write (io, '(a)') ' Matrix elements for the xc potential'
      call twocenter_write_header(io, dmax, .false.)
      close(io)
      if (ispec1_ /= ispec2_) then
        io = 3620
        write (fname2, '(a,i2.2,a,i2.2,a)') 'vxc_2c.', nzx2_, '.',             &
        &     nzx1_, '.dat'
        write (stdout, '(a)') 'Computing '//trim(fname2)//'...'
        open (unit=io, file='coutput/'//trim(fname2), status='new',            &
        &     action='write', iostat=twocenter_calc)
        if (twocenter_calc /= 0) then
          write (stderr, '(a)') '[ERROR] failed to open '//trim(fname2)
          return
        end if
        write (io, '(a)') ' Matrix elements for the xc potential'
        call twocenter_write_header(io, dmax, .true.)
        close(io)
      end if
      do issh = 1, nssh1
        io = 3710 + issh
        write (fnamel1(issh), '(a,i2.2,a,i2.2,a,i2.2,a)') 'vxc_2c_ol_', issh,  &
        &     '.', nzx1_, '.', nzx2_, '.dat'
        write (stdout, '(a)') 'Computing '//trim(fnamel1(issh))//'...'
        open (unit=io, file='coutput/'//trim(fnamel1(issh)), status='new',     &
        &     action='write', iostat=twocenter_calc)
        if (twocenter_calc /= 0) then
          write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamel1(issh))
          return
        end if
        write (io, '(a)') ' Matrix elements for the xc potential derivatives (left atom)'
        call twocenter_write_header(io, dmax, .false.)
        close(io)
        if (ispec1_ /= ispec2_) then
          io = 3720 + issh
          write (fnamer2(issh), '(a,i2.2,a,i2.2,a,i2.2,a)') 'vxc_2c_or_',      &
          &     issh, '.', nzx2_, '.', nzx1_, '.dat'
          write (stdout, '(a)') 'Computing '//trim(fnamer2(issh))//'...'
          open (unit=io, file='coutput/'//trim(fnamer2(issh)), status='new',   &
          &     action='write', iostat=twocenter_calc)
          if (twocenter_calc /= 0) then
            write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamer2(issh))
            return
          end if
          write (io, '(a)') ' Matrix elements for the xc potential derivatives (right atom)'
          call twocenter_write_header(io, dmax, .true.)
          close(io)
        end if
      end do
      do issh = 1, nssh2
        io = 3810 + issh
        write (fnamer1(issh), '(a,i2.2,a,i2.2,a,i2.2,a)') 'vxc_2c_or_', issh,  &
        &     '.', nzx1_, '.', nzx2_, '.dat'
        write (stdout, '(a)') 'Computing '//trim(fnamer1(issh))//'...'
        open (unit=io, file='coutput/'//trim(fnamer1(issh)), status='new',     &
        &     action='write', iostat=twocenter_calc)
        if (twocenter_calc /= 0) then
          write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamer1(issh))
          return
        end if
        write (io, '(a)') ' Matrix elements for the xc potential derivatives (right atom)'
        call twocenter_write_header(io, dmax, .false.)
        close(io)
        if (ispec1_ /= ispec2_) then
          io = 3820 + issh
          write (fnamel2(issh), '(a,i2.2,a,i2.2,a,i2.2,a)') 'vxc_2c_ol_',      &
          &     issh, '.', nzx2_, '.', nzx1_, '.dat'
          write (stdout, '(a)') 'Computing '//trim(fnamel2(issh))//'...'
          open (unit=io, file='coutput/'//trim(fnamel2(issh)), status='new',   &
          &     action='write', iostat=twocenter_calc)
          if (twocenter_calc /= 0) then
            write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamel2(issh))
            return
          end if
          write (io, '(a)') ' Matrix elements for the xc potential derivatives (left atom)'
          call twocenter_write_header(io, dmax, .true.)
          close(io)
        end if
      end do

      ! TODO: this is a disaster but for now is what we have
      ! THIS MUST BE DONE ON THE MASTER
      ! MAJOR RECODING NEEDS TO HAPPEN BEFORE THIS IS AVOIDED
      ndq = 2
      npts1 = ndq**nssh1
      npts2 = ndq**nssh2
      npts = npts1*npts2
      allocate(dq1(nssh1, ndq), dq2(nssh2, ndq), xmatt(nssh1+nssh2+1, npts),   &
      &        vvxc(npts, index_max_), stat=twocenter_calc)
      if (twocenter_calc /= 0) return
      call twocenter_charge_grid(dq1, dq2)

      d = 0.0_wp
      do igrid = 1, ndd_
        call utils_progress_bar(igrid, ndd_, 60, stdout)
        d = d + dr
        ix = 1
        do ix1 = 1, npts1
          do ix2 = 1, npts2
            call utils_progress_bar(ix, npts, 30, stdout)
            xmatt(1, ix) = 1.0_wp
            do issh = 1, nssh1
              ddq = dq1(issh, 1 + mod(ix1 - 1, ndq**issh)/ndq**(issh - 1))
              xnocc_in(issh, 1) = xnocc(issh, ispec1_) + ddq
              xmatt(issh+1, ix) = ddq
            end do
            do issh = 1, nssh2
              ddq = dq2(issh, 1 + mod(ix2 - 1, ndq**issh)/ndq**(issh - 1))
              xnocc_in(issh, 2) = xnocc(issh, ispec2_) + ddq
              xmatt(issh+1+nssh1, ix) = ddq
            end do
            twocenter_calc = twocenter_rho2c(d, xnocc_in)
            if (twocenter_calc /= 0) return
            do index = 1, index_max_
              n1 = nleft_(index)
              l1 = lleft_(index)
              m1 = mleft_(index)
              n2 = nright_(index)
              l2 = lright_(index)
              m2 = mright_(index)
              ! We call a very stripped version of the integral function
              ! The idea is to reduce the computation at a minimum because
              ! we will be making lots of them.
              vvxc(ix, index) = twocenter_integral_fit(n1, l1, m1, n2, l2, m2, d)
            end do ! index
            ix = ix + 1
          end do ! ix2
        end do ! ix1
        call utils_clean_progress_bar(30, stdout)
        ! End of the disaster

        ! Now we do the fit
        twocenter_calc =  math_lstsq(xmatt, vvxc, is_a_trans=.true.)
        if (twocenter_calc /= 0) then
          deallocate(xmatt, dq1, dq2, vvxc, fnamel1, fnamel2, fnamer1, fnamer2)
          write (stderr, '(a,i4)') '[ERROR] twocenterxc.f90: fit failed with info = ', &
          &     twocenter_calc
          return
        end if

        ! Write to disk
        io = 3610
        open (unit=io, file='coutput/'//trim(fname1), status='old',            &
        &     action='write', position='append', iostat=twocenter_calc)
        if (twocenter_calc /= 0) then
          write (stderr, '(a)') '[ERROR] failed to open '//trim(fname1)
          return
        end if
        write (io, '(4ES18.8)') (vvxc(1, index), index = 1, index_max_)
        close (io)
        if (ispec1_ /= ispec2_) then
          io = 3620
          open (unit=io, file='coutput/'//trim(fname2), status='old',           &
          &     action='write', position='append', iostat=twocenter_calc)
          if (twocenter_calc /= 0) then
            write (stderr, '(a)') '[ERROR] failed to open '//trim(fname2)
            return
          end if
          write (io, '(4ES18.8)') (vvxc(1, index), index = 1, index_max_)
          close (io)
        end if
        do issh = 1, nssh1
          io = 3710 + issh
          open (unit=io, file='coutput/'//trim(fnamel1(issh)), status='old',   &
          &     action='write', position='append', iostat=twocenter_calc)
          if (twocenter_calc /= 0) then
            write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamel1(issh))
            return
          end if
          write (io, '(4ES18.8)') (vvxc(issh + 1, index), index = 1, index_max_)
          close (io)
          if (ispec1_ /= ispec2_) then
            io = 3720 + issh
            open (unit=io, file='coutput/'//trim(fnamer2(issh)), status='old', &
            &     action='write', position='append', iostat=twocenter_calc)
            if (twocenter_calc /= 0) then
              write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamer2(issh))
              return
            end if
            write (io, '(4ES18.8)') (vvxc(issh + 1, index), index = 1, index_max_)
            close (io)
          end if
        end do
        do issh = 1, nssh2
          io = 3810 + issh
          open (unit=io, file='coutput/'//trim(fnamer1(issh)), status='old',   &
          &     action='write', position='append', iostat=twocenter_calc)
          if (twocenter_calc /= 0) then
            write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamer1(issh))
            return
          end if
          write (io, '(4ES18.8)') (vvxc(nssh1 + issh + 1, index), index = 1, index_max_)
          close (io)
          if (ispec1_ /= ispec2_) then
            io = 3820 + issh
            open (unit=io, file='coutput/'//trim(fnamel2(issh)), status='old', &
            &     action='write', position='append', iostat=twocenter_calc)
            if (twocenter_calc /= 0) then
              write (stderr, '(a)') '[ERROR] failed to open '//trim(fnamel2(issh))
              return
            end if
            write (io, '(4ES18.8)') (vvxc(nssh1 + issh + 1, index), index = 1, index_max_)
            close (io)
          end if
        end do
      end do ! grid
      call utils_clean_progress_bar(60, stdout)
      deallocate(dq1, dq2, vvxc, xmatt, fnamel1, fnamel2, fnamer1, fnamer2)
    end if

    write (stdout, '(a)') 'Done'
    twocenter_calc = 0
    return
  end function twocenter_calc

  real(wp) function twocenter_integral_fit (n1, l1, m1, n2, l2, m2, d)
    use math, only: math_Ylm
    include '../parameters.inc'
    include '../exchange.inc'
    include '../wavefunctions.inc'
    integer, intent(in) :: n1, l1, m1, n2, l2, m2
    real(kind=wp), intent(in) :: d
    integer :: irho, iz
    real(kind=wp) :: factor, factor2, factor4, phifactor, psi1, psi2,      &
    &                r1, r2, rho, rho2, z1, z2, z12, z22
    real(kind=wp), external :: vxc, psiofr, rescaled_psi
    real(kind=wp), allocatable :: rhomult(:), zmult(:)

    twocenter_integral_fit = 0.0_wp

    zmin = max(-rcutoff1_, d - rcutoff2_)
    zmax = min(rcutoff1_, d + rcutoff2_)
    dz = dz_
    nnz = int((zmax - zmin)/dz)
    if (iand(nnz, 1) == 0) nnz = nnz + 1
    rhomin = rhomin_
    rhomax = rhomax_
    drho = drho_
    nnrho = nnrho_

    allocate(zmult(nnz))
    factor = 0.33333333333333333333_wp*dz
    factor2 = 2.0_wp*factor
    factor4 = 2.0_wp*factor2
    zmult(1) = factor
    do iz = 2, nnz - 3, 2
      zmult(iz) = factor4
      zmult(iz+1) = factor2
    end do
    if (nnz > 1) zmult(nnz - 1) = factor2
    zmult(nnz) = factor

    allocate(rhomult(nnrho))
    factor = 0.33333333333333333333_wp*drho
    factor2 = 2.0_wp*factor
    factor4 = 2.0_wp*factor2
    rhomult(1) = factor
    do irho = 2, nnrho - 3, 2
      rhomult(irho) = factor4
      rhomult(irho+1) = factor2
    end do
    rhomult(nnrho_ - 1) = factor2
    rhomult(nnrho_) = factor

    do iz = 1, nnz
      z1 = zmin + real(iz - 1, kind=wp)*dz
      z2 = z1 - d
      z12 = z1*z1
      z22 = z2*z2
      do irho = 1, nnrho_
        rho = rhomin_ + real(irho - 1, kind=wp)*drho
        rho2 = rho*rho
        r1 = sqrt(z12 + rho2)
        if (r1 >= rcutoff1_) cycle
        r2 = sqrt(z22 + rho2)
        factor = zmult(iz)*rhomult(irho)
        psi1 = psiofr(ispec1_, n1, r1)
        psi1 = rescaled_psi(l1, m1, rho, r1, z1, psi1)
        psi2 = psiofr(ispec2_, n2, r2)
        psi2 = rescaled_psi(l2, m2, rho, r2, z2, psi2)
        twocenter_integral_fit = twocenter_integral_fit +                      &
        &                        vxc(rho, z1, iexc_, fraction_)*psi1*psi2*rho*factor
      end do
    end do
    phifactor = 0.5_wp*math_Ylm(l1, m1, 0.0_wp)*math_Ylm(l2, m2, 0.0_wp)
    twocenter_integral_fit = twocenter_integral_fit*phifactor
    deallocate(zmult, rhomult)
    return
  end function twocenter_integral_fit

  integer function twocenter_rho2c (d, xnocc_in)
    use constants, only: inv4pi
    include '../parameters.inc'
    include '../exchange.inc'
    include '../quadrature.inc'
    include '../wavefunctions.inc'
    real(kind=wp), intent(in) :: d, xnocc_in(:,:)
    integer :: irho, issh, iz
    real(kind=wp) :: dens, dzraw, r, r1, r2, rho, rho2, z1, z2, z12, z22
    real(kind=wp), external :: psiofr

    zmin = max(-rcutoff1_, d - rcutoff2_)
    zmax = min(rcutoff1_, d + rcutoff2_)
    rhomin = 0.0_wp
    rhomax = min(rcutoff1_, rcutoff2_)
    dzraw = min(drr_rho(ispec1_),drr_rho(ispec2_))
    dz = dzraw*ixcgridfactor
    if (dz > 0.05_wp) dz = 0.05_wp
    if (dz < 0.002_wp) dz = 0.002_wp
    nnz = int((zmax - zmin)/dz) + 1
    drho = dz
    nnrho = int((rhomax - rhomin)/drho) + 1
    if (nnrho > nrho_points) then
      write (stderr, '(a)') '[ERROR] twocenterxc.f90: nnrho > nrho_points'
      twocenter_rho2c = 1
      return
    end if
    if (nnz > nz_points) then
      write (stderr, '(a)') '[ERROR] twocenterxc.f90: nnz > nz_points'
      twocenter_rho2c = 1
      return
    end if

    do iz = 1, nnz
      z1 = zmin + (iz-1)*dz
      z2 = z1 - d
      z12 = z1*z1
      z22 = z2*z2
      do irho = 1, nnrho
        rho = rhomin + (irho-1)*drho
        rho2 = rho*rho
        r1 = sqrt(z12 + rho2)
        r2 = sqrt(z22 + rho2)
        dens = 0.0_wp
        do issh = 1, nsshxc(ispec1_)
          dens = dens + xnocc_in(issh,1)*psiofr(ispec1_,issh,r1)**2
        end do
        do issh = 1, nsshxc(ispec2_)
          dens = dens + xnocc_in(issh,2)*psiofr(ispec2_,issh,r2)**2
        end do
        rho2c(irho,iz) = dens*inv4pi
      end do
    end do

    twocenter_rho2c = 0
    if (iexc_ /= 4 .and. iexc_ /= 5 .and. iexc_ /= 6 .and. iexc_ /= 9 .and. iexc_ /= 10) return
    twocenter_rho2c = 1

    do iz = 1, nnz
      do irho = 2, nnrho - 1
        rhop2c(irho,iz) = (rho2c(irho+1,iz) - rho2c(irho-1,iz))/(2.0_wp*drho)
        rhopp2c(irho,iz) = (rho2c(irho+1,iz) - 2.0_wp*rho2c(irho,iz) + rho2c(irho-1,iz))/(drho**2)
      end do
      rhop2c(1,iz) = 2.0_wp*rhop2c(2,iz) - rhop2c(3,iz)
      rhop2c(nnrho,iz) = 2.0_wp*rhop2c(nnrho-1,iz) - rhop2c(nnrho-2,iz)
      rhopp2c(1,iz) = 2.0_wp*rhopp2c(2,iz) - rhopp2c(3,iz)
      rhopp2c(nnrho,iz) = 2.0_wp*rhopp2c(nnrho-1,iz) - rhopp2c(nnrho-2,iz)
    end do
    do irho = 1, nnrho
      do iz = 2, nnz - 1
        rhoz2c(irho,iz) = (rho2c(irho,iz+1) - rho2c(irho,iz-1))/(2.0d0*dz)
        rhozz2c(irho,iz) = (rho2c(irho,iz+1) - 2.0d0*rho2c(irho,iz) + rho2c(irho,iz-1))/(dz**2)
      end do
      rhoz2c(irho,1) = 2.0d0*rhoz2c(irho,2) - rhoz2c(irho,3)
      rhoz2c(irho,nnz) = 2.0d0*rhoz2c(irho,nnz-1) - rhoz2c(irho,nnz-2)
      rhozz2c(irho,1) = 2.0d0*rhozz2c(irho,2) - rhozz2c(irho,3)
      rhozz2c(irho,nnz) = 2.0d0*rhozz2c(irho,nnz-1) - rhozz2c(irho,nnz-2)
    end do
    do irho = 1, nnrho
      do iz = 2, nnz - 1
        rhopz2c(irho,iz) = (rhop2c(irho,iz+1) - rhop2c(irho,iz-1))/(2.0d0*dz)
      end do
      rhopz2c(irho,1) = 2.0d0*rhopz2c(irho,2) - rhopz2c(irho,3)
      rhopz2c(irho,nnz) = 2.0d0*rhopz2c(irho,nnz-1) - rhopz2c(irho,nnz-2)
    end do

    twocenter_rho2c = 0
    return
  end function twocenter_rho2c
end module twocenterxc

