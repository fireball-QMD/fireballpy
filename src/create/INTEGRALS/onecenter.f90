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


! onecenter.f90
! Program Description
! ==============================================================================
!       This module calculates the one-center integrals.
! ==============================================================================
! Original code from Juergen Fritsch
!
! Code rewritten by:
! Cyra Roldan Pinero
! ==============================================================================
module onecenter
  use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit, stdout => output_unit
  use :: constants, only: sqinv4pi, path_len, tolerance
  use :: indices, only: indices_onecenter_set
  use :: utils, only: utils_open
  use :: xc, only: xc_calc, xc_isgga
  use :: wavefunctions, only: wf_atoms, wf_dens
  implicit none
  private
  public :: onecenter_calc

  integer, parameter :: ONECENTER_NUM_INTERACTIONS = 3
  integer, parameter, public :: ONECENTER_XC        = ishft(1, 0) ! 2^0
  integer, parameter, public :: ONECENTER_GOVERLAP1 = ishft(1, 1) ! 2^1
  integer, parameter, public :: ONECENTER_GOVERLAP2 = ishft(1, 2) ! 2^2

contains

  subroutine onecenter_get_fname(interaction, ispec, fname)
    integer, intent(in) :: interaction, ispec
    character(path_len), intent(out) :: fname
    character(32) :: root
    character(2) :: auxz
    select case (interaction)
    case (ONECENTER_XC)
      root = "xc"
    case (ONECENTER_GOVERLAP1)
      root = "goverlapf1"
    case (ONECENTER_GOVERLAP2)
      root = "goverlapf2"
    case default
      error stop
    end select
    write (auxz,"(i2.2)") wf_atoms(ispec)%get_nz()
    fname = trim(root)//"."//auxz//".dat"
  end subroutine onecenter_get_fname

  pure integer function onecenter_get_nints(interaction, ispec)
    integer, intent(in) :: interaction, ispec
    integer :: nsh
    nsh = wf_atoms(ispec)%get_nshells()
    select case (interaction)
    case (ONECENTER_XC)
      onecenter_get_nints = 2 + 2*nsh
    case default
      onecenter_get_nints = 1
    end select
  end function onecenter_get_nints

  pure subroutine onecenter_get_interactions(interactions, onecenter_interactions)
    integer, intent(in) :: interactions
    logical, intent(out) :: onecenter_interactions(ONECENTER_NUM_INTERACTIONS)
    integer :: i, p
    onecenter_interactions = .false.
    p = ishft(1, ONECENTER_NUM_INTERACTIONS - 1)
    do i = ONECENTER_NUM_INTERACTIONS, 1, -1
      if (iand(interactions, p) == p) then
        onecenter_interactions(i) = .true.
      end if
      p = ishft(p, -1)
    end do
  end subroutine onecenter_get_interactions

  subroutine onecenter_calc_interaction(interaction, ispec, answer)
    integer, intent(in) :: interaction, ispec
    real(dp), allocatable, intent(out) :: answer(:,:)
    integer :: irho, ish, jssh, iint, index, nsh, nints, index_max, nrho
    real(kind=dp) :: ir, rcut, factor, drho, rho, tmp, dtmp, ddtmp, psi1, psi2, dens, lapl, &
      &              exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross
    logical :: onecenter_interactions(ONECENTER_NUM_INTERACTIONS)
    integer, allocatable :: s1(:), s2(:), l12(:)
    real(kind=dp), allocatable :: fofr(:), grad(:), hess(:,:)

    ! Retrieve basic info
    nsh = wf_atoms(ispec)%get_nshells()
    drho = wf_atoms(ispec)%get_dr()
    rcut = wf_atoms(ispec)%get_rcut()
    nrho = rcut/drho + 1
    if (iand(nrho, 1) == 0) then
      nrho = nrho + 1
      drho = rcut/real(nrho - 1, kind=dp)
    end if

    ! Set dimensions
    nints = onecenter_get_nints(interaction, ispec)
    call indices_onecenter_set([(wf_atoms(ispec)%get_angular_momentum(ish), ish = 1, nsh)], index_max, s1, s2, l12)
    if (allocated(answer)) deallocate(answer)
    allocate(fofr(nints), answer(nints, index_max))
    answer = 0.0_dp

    do irho = 2, nrho ! We can ignore 0 because of the factor term
      rho = real(irho - 1, kind=dp)*drho
      factor = 0.66666666666666666667_dp*drho
      if (iand(irho, 1) == 0) factor = 2.0_dp*factor
      if (irho == nrho) factor = 0.5_dp*factor

      ! The integral in all its glory
      select case (interaction)
      case (ONECENTER_XC)
        if (xc_isgga()) then
          call wf_dens(ispec, rho, dens, lapl, grad, hess)
          call xc_calc(dens, lapl, grad, hess, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma, dvxclapl, dvxccross)
          if (rho > tolerance) ir = 1.0_dp/rho
        else
          call wf_dens(ispec, rho, dens)
          call xc_calc(dens, exc, vxc, dexcrho, dvxcrho)
        end if
        fofr(1) = vxc
        fofr(nsh + 2) = exc
        do iint = 1, nsh
          tmp = sqinv4pi*wf_atoms(ispec)%get_psi(iint, rho)
          fofr(iint + 1) = dvxcrho*tmp*tmp
          fofr(iint + 2 + nsh) = dexcrho*tmp*tmp
          if (xc_isgga() .and. rho > tolerance) then
            dtmp = sqinv4pi*wf_atoms(ispec)%get_psi(iint, rho, order=1)
            ddtmp = sqinv4pi*wf_atoms(ispec)%get_psi(iint, rho, order=2)
            fofr(iint + 1) = fofr(iint + 1) + dvxcsigma*4.0_dp*tmp*dtmp*grad(1) + &
              &                               dvxclapl*2.0_dp*(dtmp*dtmp + tmp*ddtmp + 2.0_dp*ir*tmp*dtmp) + &
              &                               dvxccross*2.0_dp*grad(1)*(grad(1)*(dtmp*dtmp + tmp*ddtmp) + 2.0_dp*tmp*dtmp*hess(1,1))
            fofr(iint + 2 + nsh) = fofr(iint + 2 + nsh) + dexcsigma*4.0_dp*tmp*dtmp*grad(1)
          end if
        end do
      case default
        fofr(1) = 1.0_dp
      end select

      do index = 1, index_max
        ! Second wavefunction changes
        psi1 = wf_atoms(ispec)%get_psi(s1(index), rho)
        select case (interaction)
        case (ONECENTER_GOVERLAP1)
          psi2 = wf_atoms(ispec)%get_psi(s2(index), rho, order=1)
        case default
          psi2 = wf_atoms(ispec)%get_psi(s2(index), rho)
        end select

        ! Multiply by the psi's
        select case (interaction)
        case (ONECENTER_GOVERLAP2)
          answer(:, index) = answer(:, index) + fofr(:)*psi1*psi2*factor*rho
        case default
          answer(:, index) = answer(:, index) + fofr(:)*psi1*psi2*factor*rho*rho
        end select
      end do
    end do
  end subroutine onecenter_calc_interaction

  subroutine onecenter_write_interaction(interaction, ispec, fname, answer)
    integer, intent(in) :: interaction, ispec
    character(path_len), intent(in) :: fname
    real(dp), intent(in) :: answer(:,:)
    integer :: ish, jssh, ix, index, index_max, nsh, nzx, io, nints
    integer, allocatable :: s1(:), s2(:), l12(:)
    real(dp) :: rcut
    real(dp), allocatable :: matrix(:,:)
    character(64), allocatable :: names(:)

    ! Retrieve basic info
    nints = onecenter_get_nints(interaction, ispec)
    index_max = size(answer, dim=2)
    nsh = wf_atoms(ispec)%get_nshells()
    rcut = wf_atoms(ispec)%get_rcut()
    nzx = wf_atoms(ispec)%get_nz()

    io = utils_open(fname, "w")
    write (io, "(14x,i2,44x,'! Atomic number')") nzx
    write (io, "(2x,f14.6,44x,'! Cutoff radius')") rcut
    if (interaction == ONECENTER_XC) then
      write (io, "(8x,i8,36x,'! Number of shells')") nsh
      write (io, "(1000ES16.8)", advance="no") (wf_atoms(ispec)%get_ref_charge(ish), ish = 1, nsh)
      write (io, "(a)") " ! Reference charges"
    end if

    ! Create the matrix to output in matrix format
    call indices_onecenter_set([(wf_atoms(ispec)%get_angular_momentum(ish), ish = 1, nsh)], &
      &                        index_max, s1, s2, l12, names)
    write (io, "('!')", advance="no")
    if (index_max == 0) then
      close (io)
      return
    end if
    do index = 1, index_max
      write (io, "(4x,a)", advance="no") trim(names(index))
    end do
    write (io, "(a)") ""
    allocate(matrix(nsh, nsh))
    do ix = 1, nints
      matrix = 0.0_dp
      do index = 1, index_max
        matrix(s1(index), s2(index)) = answer(ix, index)
      end do
      do ish = 1, nsh
        write (io, "(1000ES16.8)") (matrix(ish, jssh), jssh = 1, nsh)
      end do
      if (ix < nints) write (io, "(a)") ""
    end do
    close (io)
  end subroutine onecenter_write_interaction

  subroutine onecenter_calc(interactions, ispec)
    integer, intent(in) :: interactions, ispec
    integer :: interaction, int_id, iint, nints
    logical :: exists
    logical :: onecenter_interactions(ONECENTER_NUM_INTERACTIONS)
    character(path_len) :: fname
    real(dp), allocatable :: answer(:,:)

    call onecenter_get_interactions(interactions, onecenter_interactions)
    do int_id = 1, ONECENTER_NUM_INTERACTIONS
      if (.not. onecenter_interactions(int_id)) cycle
      interaction = ishft(1, int_id - 1)
      call onecenter_get_fname(interaction, ispec, fname)
      inquire (file=fname, exist=exists)
      if (exists) cycle
      write (stdout, "(2x,a)") "Computing "//trim(fname)//"..."
      call onecenter_calc_interaction(interaction, ispec, answer)
      call onecenter_write_interaction(interaction, ispec, fname, answer)
      write (stdout, "(a)", advance="no") achar(8)//achar(13)
      write (stdout, "(2x,a)") "Computing "//trim(fname)//"... Done!"
    end do ! int_id
  end subroutine onecenter_calc

end module onecenter
