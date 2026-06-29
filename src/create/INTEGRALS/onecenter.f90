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
  use :: constants, only: inv4pi, path_len
  use :: indices, only: indices_onecenter_set
  use :: utils, only: utils_open
  use :: xc, only: xc_calc, xc_isgga
  use :: wavefunctions, only: wf_atoms, wf_dens
  implicit none
  private
  public :: onecenter_calc

  integer, parameter, public :: ONECENTER_NPOINTS = 1025

  integer, parameter :: ONECENTER_NUM_INTERACTIONS = 3
  integer, parameter, public :: ONECENTER_XC        = ishft(1, 0) ! 2^0
  integer, parameter, public :: ONECENTER_GOVERLAP1 = ishft(1, 1) ! 2^1
  integer, parameter, public :: ONECENTER_GOVERLAP2 = ishft(1, 2) ! 2^2
  character(10), parameter :: ONECENTER_ROOT(ONECENTER_NUM_INTERACTIONS) = [ &
    &  "xc        ", &
    &  "goverlapf1", &
    &  "goverlapf2" &
    &  ]

contains

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

  pure integer function onecenter_get_nints(int_id, ispec)
    integer, intent(in) :: int_id, ispec
    integer :: interaction, nssh
    interaction = ishft(1, int_id - 1)
    nssh = wf_atoms(ispec)%get_nshells()
    select case (interaction)
    case (ONECENTER_XC)
      onecenter_get_nints = 2 + 2*nssh
    case default
      onecenter_get_nints = 1
    end select
  end function onecenter_get_nints

  subroutine onecenter_calc_interaction(int_id, ispec, answer)
    integer, intent(in) :: int_id, ispec
    real(dp), allocatable, intent(out) :: answer(:,:)
    integer :: interaction, irho, issh, jssh, isorp, index, nssh, nints, index_max
    real(kind=dp) :: rcut, factor, drho, rho, tmp, psi1, psi2
    logical :: onecenter_interactions(ONECENTER_NUM_INTERACTIONS)
    integer, allocatable :: s1(:), s2(:)
    real(kind=dp), allocatable :: exc(:), vxc(:), dexcrho(:), dexcsigma(:), dvxcrho(:), dvxcsigma(:), fofr(:), &
      &                           dens(:), ddens(:,:), dddens(:,:,:)

    interaction = ishft(1, int_id - 1)

    ! Retrieve basic info
    nssh = wf_atoms(ispec)%get_nshells()
    rcut = wf_atoms(ispec)%get_rcut()

    ! Set dimensions
    nints = onecenter_get_nints(int_id, ispec)
    call indices_onecenter_set([(wf_atoms(ispec)%get_angular_momentum(issh), issh = 1, nssh)], &
      &                        index_max, s1, s2)
    allocate (fofr(nints), answer(nints, index_max))
    answer = 0.0_dp

    ! Allocate and compute for XC
    if (interaction == ONECENTER_XC) then
      allocate (dens(ONECENTER_NPOINTS), exc(ONECENTER_NPOINTS), vxc(ONECENTER_NPOINTS), &
        &       dexcrho(ONECENTER_NPOINTS), dvxcrho(ONECENTER_NPOINTS))
      if (xc_isgga()) then
        allocate (ddens(ONECENTER_NPOINTS, 1), dddens(ONECENTER_NPOINTS, 1, 1), &
        &         dexcsigma(ONECENTER_NPOINTS), dvxcsigma(ONECENTER_NPOINTS))
        call wf_dens(ispec, dens, ddens, dddens)
        call xc_calc(dens, ddens, dddens, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma)
      else
        call wf_dens(ispec, dens)
        call xc_calc(dens, exc, vxc, dexcrho, dvxcrho)
      end if
    end if

    drho = rcut/real(ONECENTER_NPOINTS - 1, kind=dp)
    do irho = 2, ONECENTER_NPOINTS  ! Avoid 0 problems
      rho = real(irho - 1, kind=dp)*drho
      factor = 0.66666666666666666667_dp*drho
      if (iand(irho, 1) == 0) factor = 2.0_dp*factor
      if (irho == 1 .or. irho == ONECENTER_NPOINTS) factor = 0.5_dp*factor

      ! The integral in all its glory
      select case (interaction)
      case (ONECENTER_XC)
        fofr(1) = vxc(irho)
        do isorp = 1, nssh
          tmp = wf_atoms(ispec)%get_psi(isorp, rho)
          fofr(isorp + 1) = inv4pi*tmp*tmp*dvxcrho(irho)
          if (xc_isgga()) then
            fofr(isorp + 1) = fofr(isorp + 1) + 4.0_dp*inv4pi*tmp*dvxcsigma(irho) * &
              &               wf_atoms(ispec)%get_psi(isorp, rho, order=1)*ddens(irho, 1)
          end if
        end do
        fofr(nssh + 2) = exc(irho)
        do isorp = 1, nssh
          tmp = wf_atoms(ispec)%get_psi(isorp, rho)
          fofr(isorp + 2 + nssh) = inv4pi*tmp*tmp*dexcrho(irho)
          if (xc_isgga()) then
            fofr(isorp + 2 + nssh) = fofr(isorp + 2 + nssh) + 4.0_dp*inv4pi*tmp*dexcsigma(irho) * &
              &                      wf_atoms(ispec)%get_psi(isorp, rho, order=1)*ddens(irho, 1)
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
    deallocate (s1, s2, fofr)

    ! Don't forget to clean XC
    if (interaction == ONECENTER_XC) then
      deallocate (dens, exc, vxc, dexcrho, dvxcrho)
      if (xc_isgga()) deallocate (ddens, dddens, dexcsigma, dvxcsigma)
    end if
  end subroutine onecenter_calc_interaction

  pure subroutine onecenter_get_fname(int_id, ispec, fname)
    integer, intent(in) :: int_id, ispec
    character(path_len), intent(out) :: fname
    integer :: nzx
    character(2) :: auxz

    ! Retrieve basic info
    nzx = wf_atoms(ispec)%get_nz()
    write (auxz,"(i2.2)") nzx
    fname = trim(ONECENTER_ROOT(int_id))//"."//auxz//".dat"
  end subroutine onecenter_get_fname

  subroutine onecenter_write_interaction(int_id, ispec, fname, answer)
    integer, intent(in) :: int_id, ispec
    character(path_len), intent(in) :: fname
    real(dp), intent(in) :: answer(:,:)
    integer :: issh, jssh, ix, index, index_max, nssh, nzx, interaction, io, nints
    integer, allocatable :: s1(:), s2(:)
    real(dp) :: rcut
    real(dp), allocatable :: matrix(:,:)

    ! Retrieve basic info
    nints = onecenter_get_nints(int_id, ispec)
    index_max = size(answer, dim=2)
    nssh = wf_atoms(ispec)%get_nshells()
    rcut = wf_atoms(ispec)%get_rcut()
    nzx = wf_atoms(ispec)%get_nz()
    interaction = ishft(1, int_id - 1)

    io = utils_open(fname, "w")
    write (io, "(14x,i2,44x,'! Atomic number')") nzx
    write (io, "(2x,f14.6,44x,'! Cutoff radius')") rcut
    if (interaction == ONECENTER_XC) then
      write (io, "(8x,i8,36x,'! Number of shells')") nssh
      write (io, "(1000ES16.8)") (wf_atoms(ispec)%get_ref_charge(issh), issh = 1, nssh)
    end if

    ! Create the matrix to output in matrix format
    call indices_onecenter_set([(wf_atoms(ispec)%get_angular_momentum(issh), issh = 1, nssh)], &
      &                        index_max, s1, s2)
    allocate (matrix(nssh, nssh))
    write (io, "(a)") ""
    do ix = 1, nints
      matrix = 0.0_dp
      do index = 1, index_max
        matrix(s1(index), s2(index)) = answer(ix, index)
      end do
      do issh = 1, nssh
        write (io, "(1000ES16.8)") (matrix(issh, jssh), jssh = 1, nssh)
      end do
      if (ix < nints) write (io, "(a)") ""
    end do
    close (io)
    deallocate (s1, s2, matrix)
  end subroutine onecenter_write_interaction

  subroutine onecenter_calc(interactions, ispec)
    integer, intent(in) :: interactions, ispec
    integer :: interaction, int_id, isorp, nints
    logical :: exists
    logical :: onecenter_interactions(ONECENTER_NUM_INTERACTIONS)
    character(path_len) :: fname
    real(dp), allocatable :: answer(:,:)

    if (iand(ONECENTER_NPOINTS, 1) == 0) then
      write (stderr, "('[ERROR] ',(a))") "Number of integration points must be odd"
      error stop 1
    end if
    if (ONECENTER_NPOINTS == 1) then
      write (stderr, "('[ERROR] ',(a))") "Number of integration points must be > 1"
      error stop 1
    end if

    call onecenter_get_interactions(interactions, onecenter_interactions)
    do int_id = 1, ONECENTER_NUM_INTERACTIONS
      if (.not. onecenter_interactions(int_id)) cycle
      interaction = ishft(1, int_id - 1)
      call onecenter_get_fname(int_id, ispec, fname)
      inquire (file=fname, exist=exists)
      if (exists) cycle
      write (stdout, "(2x,a)") "Computing "//trim(fname)//"..."
      call onecenter_calc_interaction(int_id, ispec, answer)
      call onecenter_write_interaction(int_id, ispec, fname, answer)
      write (stdout, "(a)", advance="no") achar(8)//achar(13)
      write (stdout, "(2x,a)") "Computing "//trim(fname)//"... Done!"
    end do ! int_id
    deallocate (answer)
  end subroutine onecenter_calc

end module onecenter
