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
! Brigham Young University - Hao Wang Lawrence
! Livermore National Laboratory - Kurt Glaesemann
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
!       This module calculates the two-center integrals.
! ==============================================================================
! Code written by:
! Cyra Roldan Pinero
! ==============================================================================
module twocenter
  use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit, stdout => output_unit
  use :: constants, only: inv4pi, sqinv4pi, twopi, invsq2, tolerance, path_len, fname_len
  use :: utils, only: utils_open
  use :: indices, only: indices_twocenter_set, INDICES_TWOCENTER, INDICES_TWOCENTER_DIPX, &
    &                   INDICES_TWOCENTER_DIPY, INDICES_TWOCENTER_COULOMB
  use :: xc, only: xc_calc, xc_isgga
  use :: wavefunctions, only: wf_atoms, wf_dens
  use :: potentials, only: pot_atoms
  use :: pseudopotentials, only: pp_atoms
  implicit none
  private
  public :: twocenter_calc

  integer, parameter, public :: TWOCENTER_NPOINTS_D = 128
  integer, parameter, public :: TWOCENTER_NPOINTS_Z = 129
  integer, parameter, public :: TWOCENTER_NPOINTS_RHO = 129

  integer, parameter :: TWOCENTER_NUM_INTERACTIONS = 14
  integer, parameter, public :: TWOCENTER_DENS_ATOM      = ishft(1, 0)
  integer, parameter, public :: TWOCENTER_DENS_ONTOP     = ishft(1, 1)
  integer, parameter, public :: TWOCENTER_OVERLAP        = ishft(1, 2)
  integer, parameter, public :: TWOCENTER_VNA_ATOM       = ishft(1, 3)
  integer, parameter, public :: TWOCENTER_VNA_ONTOP      = ishft(1, 4)
  integer, parameter, public :: TWOCENTER_VPP            = ishft(1, 5)
  integer, parameter, public :: TWOCENTER_VXC            = ishft(1, 6)
  integer, parameter, public :: TWOCENTER_DIP_Z          = ishft(1, 7)
  integer, parameter, public :: TWOCENTER_DIP_X          = ishft(1, 8)
  integer, parameter, public :: TWOCENTER_DIP_Y          = ishft(1, 9)
  integer, parameter, public :: TWOCENTER_COULOMB        = ishft(1, 10)
  integer, parameter, public :: TWOCENTER_DENS_ATOM_SPH  = ishft(1, 11)
  integer, parameter, public :: TWOCENTER_DENS_ONTOP_SPH = ishft(1, 12)
  integer, parameter, public :: TWOCENTER_OVERLAP_SPH    = ishft(1, 13)
  character(10), parameter :: TWOCENTER_ROOT(TWOCENTER_NUM_INTERACTIONS) = [ &
    &  "den_atom  ", &
    &  "den_ontop ", &
    &  "overlap   ", &
    &  "vna_atom  ", &
    &  "vna_ontop ", &
    &  "vnl       ", &
    &  "vxc       ", &
    &  "dipole_z  ", &
    &  "dipole_x  ", &
    &  "dipole_y  ", &
    &  "coulomb   ", &
    &  "denS_atom ", &
    &  "denS_ontop", &
    &  "overlapS  " &
    &  ]

contains

  pure subroutine twocenter_get_interactions(interactions, twocenter_interactions)
    integer, intent(in) :: interactions
    logical, intent(out) :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    integer :: i, p
    twocenter_interactions = .false.
    p = ishft(1, TWOCENTER_NUM_INTERACTIONS - 1)
    do i = TWOCENTER_NUM_INTERACTIONS, 1, -1
      if (iand(interactions, p) == p) then
        twocenter_interactions(i) = .true.
      end if
      p = ishft(p, -1)
    end do
  end subroutine twocenter_get_interactions

  pure integer function twocenter_get_nints(int_id, ispec, jspec)
    integer, intent(in) :: int_id, ispec, jspec
    integer :: interaction, nsshi, nsshj
    interaction = ishft(1, int_id - 1)
    nsshi = wf_atoms(ispec)%get_nshells()
    nsshj = wf_atoms(jspec)%get_nshells()
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_DENS_ATOM_SPH)
      twocenter_get_nints = nsshj
    case (TWOCENTER_VNA_ATOM)
      twocenter_get_nints = nsshj + 1
    case (TWOCENTER_DENS_ONTOP, TWOCENTER_DENS_ONTOP_SPH)
      twocenter_get_nints = nsshi + nsshj
    case (TWOCENTER_VXC)
      twocenter_get_nints = nsshi + nsshj + 1
    case (TWOCENTER_VNA_ONTOP)
      twocenter_get_nints = nsshi + nsshj + 2
    case default
      twocenter_get_nints = 1
    end select
  end function twocenter_get_nints

  pure integer function twocenter_get_index_type(int_id)
    integer, intent(in) :: int_id
    integer :: interaction
    interaction = ishft(1, int_id - 1)
    select case (interaction)
    case (TWOCENTER_DIP_X)
      twocenter_get_index_type = INDICES_TWOCENTER_DIPX
    case (TWOCENTER_DIP_Y)
      twocenter_get_index_type = INDICES_TWOCENTER_DIPY
    case (TWOCENTER_COULOMB)
      twocenter_get_index_type = INDICES_TWOCENTER_COULOMB
    case default
      twocenter_get_index_type = INDICES_TWOCENTER
    end select
  end function twocenter_get_index_type

  subroutine twocenter_calc_interaction(int_id, ispec, jspec, index_max, answer, names)
    integer, intent(in) :: int_id, ispec, jspec
    integer, intent(out) :: index_max
    real(dp), allocatable, intent(out) :: answer(:,:,:)
    character(64), allocatable, intent(out) :: names(:)
    integer :: spec1, spec2, nssh1, nssh2, nsshi, nsshj, nints, &
      &        interaction, igrid, iz, irho, isorp, index
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    real(dp) :: d, dd, rcuti, rcutj, zmin, zmax, dz, drho, zi, zj, zi2, zj2, ri, rj, z1, z2, r1, r2, &
      &         zmult, factor, psi1, psi2, cyl1, cyl2, nzxi, nzxj, rcut1, rcut2, tmp, rho, rho2
    integer, allocatable :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:), ls1(:), ls2(:)
    real(dp), allocatable :: rhomult(:), fofr(:)
    real(dp), allocatable :: dens(:), exc(:), vxc(:), dexcrho(:), dvxcrho(:), &
      &                      dexcsigma(:), dvxcsigma(:), ddens(:,:), dddens(:,:,:)

    interaction = ishft(1, int_id - 1)

    ! Retrieve basic info
    nsshi = wf_atoms(ispec)%get_nshells()
    nsshj = wf_atoms(jspec)%get_nshells()
    rcuti = wf_atoms(ispec)%get_rcut()
    rcutj = wf_atoms(jspec)%get_rcut()

    ! First atom is always normal (save for spherical), the second adapts
    spec1 = ispec
    nssh1 = wf_atoms(spec1)%get_nshells()
    rcut1 = rcuti
    ls1 = [(wf_atoms(spec1)%get_angular_momentum(index), index = 1, nssh1)]
    select case (interaction)  ! spec2 (atom)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNA_ATOM, TWOCENTER_DENS_ATOM_SPH)
      spec2 = ispec
    case default
      spec2 = jspec
    end select
    select case (interaction)  ! nssh2, rcut2, ls2 (PP and long range)
    case (TWOCENTER_DENS_ATOM_SPH, TWOCENTER_DENS_ONTOP_SPH, TWOCENTER_OVERLAP_SPH)
      nssh2 = wf_atoms(spec2)%get_nshells()
      rcut2 = wf_atoms(spec2)%get_rcut()
      ls1 = 0
      ls2 = [(0, index = 1, nssh2)]
    case (TWOCENTER_VPP)
      nssh2 = pp_atoms(spec2)%get_nshells()
      rcut2 = pp_atoms(spec2)%get_rcut()
      ls2 = [(pp_atoms(spec2)%get_angular_momentum(index), index = 1, nssh2)]
    case (TWOCENTER_COULOMB)
      nssh2 = wf_atoms(spec2)%get_nshells()
      rcut2 = 1.0e10_dp
      ls2 = [(wf_atoms(spec2)%get_angular_momentum(index), index = 1, nssh2)]
    case default
      nssh2 = wf_atoms(spec2)%get_nshells()
      rcut2 = rcutj
      ls2 = [(wf_atoms(spec2)%get_angular_momentum(index), index = 1, nssh2)]
    end select

    ! Set dimensions
    nints = twocenter_get_nints(int_id, ispec, jspec)
    call indices_twocenter_set(twocenter_get_index_type(int_id), ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    deallocate (ls1, ls2)
    if (index_max == 0) return
    if (allocated(answer)) deallocate (answer)
    allocate (fofr(nints), answer(nints, index_max, TWOCENTER_NPOINTS_D))
    answer = 0.0_dp

    ! Allocate for XC
    if (interaction == TWOCENTER_VXC) then
      allocate (dens(TWOCENTER_NPOINTS_RHO), exc(TWOCENTER_NPOINTS_RHO), vxc(TWOCENTER_NPOINTS_RHO), &
        &       dexcrho(TWOCENTER_NPOINTS_RHO), dvxcrho(TWOCENTER_NPOINTS_RHO))
      if (xc_isgga()) allocate (ddens(TWOCENTER_NPOINTS_RHO, 2), dddens(TWOCENTER_NPOINTS_RHO, 2, 2), &
          &                     dexcsigma(TWOCENTER_NPOINTS_RHO), dvxcsigma(TWOCENTER_NPOINTS_RHO))
    end if

    ! Prepare integral
    dd = (rcut1 + rcut2)/real(TWOCENTER_NPOINTS_D - 1, kind=dp)
    drho = min(rcut1, rcut2)/real(TWOCENTER_NPOINTS_RHO - 1, kind=dp)
    allocate(rhomult(TWOCENTER_NPOINTS_RHO))
    do irho = 1, TWOCENTER_NPOINTS_RHO
      rhomult(irho) = 0.66666666666666666667_dp*drho
      if (irho == 1 .or. irho == TWOCENTER_NPOINTS_RHO) then
        rhomult(irho) = 0.5_dp*rhomult(irho)
      else if (iand(irho, 1) == 0) then
        rhomult(irho) = 2.0_dp*rhomult(irho)
      end if
    end do

    do igrid = 1, TWOCENTER_NPOINTS_D
      d = real(igrid - 1, kind=dp)*dd
      zmin = max(-rcut1, d - rcut2)
      zmax = min(rcut1, d + rcut2)
      dz = (zmax - zmin)/real(TWOCENTER_NPOINTS_Z - 1, kind=dp)

      do iz = 1, TWOCENTER_NPOINTS_Z
        zi = zmin + real(iz - 1, kind=dp)*dz
        zj = zi - d
        zi2 = zi*zi
        zj2 = zj*zj
        zmult = 0.66666666666666666667_dp*dz
        if (iz == 1 .or. iz == TWOCENTER_NPOINTS_Z) then
          zmult = 0.5_dp*zmult
        else if (iand(iz, 1) == 0) then
          zmult = 2.0_dp*zmult
        end if

        ! XC is computed at once for all rho values
        if (interaction == TWOCENTER_VXC) then
          if (xc_isgga()) then
            call wf_dens(spec1, spec2, z1, z2, dens, ddens, dddens)
            call xc_calc(dens, ddens, dddens, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma)
          else
            call wf_dens(spec1, spec2, z1, z2, dens)
            call xc_calc(dens, exc, vxc, dexcrho, dvxcrho)
          end if
        end if

        do irho = 2, TWOCENTER_NPOINTS_RHO  ! Avoid 0 problems
          rho = real(irho - 1, kind=dp)*drho
          rho2 = rho*rho
          ri = sqrt(zi2 + rho2)
          if (ri >= rcut1) cycle
          rj = sqrt(zj2 + rho2)
          if (rj >= rcut2) cycle
          factor = zmult*rhomult(irho)*twopi

          ! Atom case (again only second atom adapts)
          z1 = zi
          r1 = ri
          select case (interaction)
          case (TWOCENTER_DENS_ATOM, TWOCENTER_VNA_ATOM, TWOCENTER_DENS_ATOM_SPH)
            z2 = zi
            r2 = ri
          case default
            z2 = zj
            r2 = rj
          end select

          ! The integral in all its glory
          select case (interaction)
          case (TWOCENTER_DENS_ATOM, TWOCENTER_DENS_ATOM_SPH)
            do isorp = 1, nsshj
              fofr(isorp) = inv4pi*wf_atoms(spec2)%get_psi(isorp, rj)**2
            end do
          case (TWOCENTER_DENS_ONTOP, TWOCENTER_DENS_ONTOP_SPH)
            do isorp = 1, nsshi
              fofr(isorp) = inv4pi*wf_atoms(spec1)%get_psi(isorp, ri)**2
            end do
            do isorp = 1, nsshj
              fofr(isorp + nsshi) = inv4pi*wf_atoms(spec2)%get_psi(isorp, rj)**2
            end do
          case (TWOCENTER_VNA_ATOM)
            do isorp = 0, nsshj
              fofr(isorp + 1) = pot_atoms(spec2)%get_vnn(isorp, rj)**2
            end do
          case (TWOCENTER_VNA_ONTOP)
            do isorp = 0, nsshi
              fofr(isorp + 1) = pot_atoms(spec1)%get_vnn(isorp, ri)**2
            end do
            do isorp = 0, nsshj
              fofr(isorp + nsshi + 2) = pot_atoms(spec2)%get_vnn(isorp, rj)**2
            end do
          case (TWOCENTER_VXC)
            fofr(1) = vxc(irho)
            do isorp = 1, nsshi
              tmp = wf_atoms(spec1)%get_psi(isorp, r1)
              fofr(isorp + 1) = inv4pi*tmp*tmp*dvxcrho(irho)
              if (xc_isgga() .and. r1 > tolerance) then
                fofr(isorp + 1) = fofr(isorp + 1) + dvxcsigma(irho)*4.0_dp*inv4pi * &
                  &               tmp*wf_atoms(spec1)%get_psi(isorp, r1, order=1) * &
                  &               (ddens(irho, 1)*rho + ddens(irho, 2)*z1)/r1
              end if
            end do
            do isorp = 1, nsshj
              tmp = wf_atoms(spec2)%get_psi(isorp, r2)
              fofr(isorp + 1 + nsshi) = inv4pi*tmp*tmp*dvxcrho(irho)
              if (xc_isgga() .and. r2 > tolerance) then
                fofr(isorp + 1 + nsshi) = fofr(isorp + 1 + nsshi) + dvxcsigma(irho)*4.0_dp*inv4pi * &
                  &                       tmp*wf_atoms(spec2)%get_psi(isorp, r2, order=1) * &
                  &                       (ddens(irho, 1)*rho + ddens(irho, 2)*z2)/r2
              end if
            end do
          case (TWOCENTER_DIP_Z)
            fofr(1) = z1 - 0.5_dp*d
          case (TWOCENTER_DIP_X, TWOCENTER_DIP_Y)
            fofr(1) = rho
          case default
            fofr(1) = 1.0_dp
          end select

          do index = 1, index_max
            ! Second atom wavefunction changes
            psi1 = wf_atoms(spec1)%get_psi(s1(index), r1)
            select case (interaction)
            case (TWOCENTER_VPP)
              psi2 = pp_atoms(spec2)%get_vpp(s2(index), r2)
            case (TWOCENTER_COULOMB)
              psi2 = pot_atoms(spec2)%get_vnn(s2(index), r2)
            case default
              psi2 = wf_atoms(spec2)%get_psi(s2(index), r2)
            end select

            ! Get the cylindrical factors for the wavefunctions (not Coulomb)
            select case (interaction)
            case (TWOCENTER_COULOMB)
              cyl1 = 1.0_dp
              cyl2 = 1.0_dp
            case default
              cyl1 = cylindrical_harmonic(l1(index), m1(index), rho, z1, r1)
              cyl2 = cylindrical_harmonic(l2(index), m2(index), rho, z2, r2)
            end select

            ! Some signs and factors need to be adjusted for dipX and dipY
            if (interaction == TWOCENTER_DIP_X .or. interaction == TWOCENTER_DIP_Y) then
              if (iand(m1(index), 1) == 0) cyl1 = cyl1*invsq2
              if (iand(m2(index), 1) == 0) cyl2 = cyl2*invsq2
              if (interaction == TWOCENTER_DIP_Y) then
                if (l1(index) == 2 .and. m1(index) == 2 .or. &
                  & l1(index) == 3 .and. (m1(index) == -3 .or. m1(index) == 2 .or. m1(index) == 3)) cyl1 = -cyl1
                if (l2(index) == 2 .and. m2(index) == 2 .or. &
                  & l2(index) == 3 .and. (m2(index) == -3 .or. m2(index) == 2 .or. m2(index) == 3)) cyl2 = -cyl2
              end if
            end if

            ! Multiply by the psi's
            select case (interaction)
            case (TWOCENTER_DENS_ATOM_SPH, TWOCENTER_DENS_ONTOP_SPH, TWOCENTER_OVERLAP_SPH)
              answer(:, index, igrid) = answer(:, index, igrid) + fofr(:)*abs(psi1*psi2)*factor*rho
            case default
              answer(:, index, igrid) = answer(:, index, igrid) + fofr(:)*psi1*psi2*cyl1*cyl2*factor*rho
            end select
          end do ! index
        end do ! irho
      end do ! iz
    end do ! igrid
    deallocate (s1, s2, l1, l2, m1, m2, fofr)

    ! Don't forget to clean XC
    if (interaction == TWOCENTER_VXC) then
      deallocate (dens, exc, vxc, dexcrho, dvxcrho)
      if (xc_isgga()) deallocate (ddens, dddens, dexcsigma, dvxcsigma)
    end if
  end subroutine twocenter_calc_interaction

  subroutine twocenter_get_fnames(int_id, ispec, jspec, fnames)
    integer, intent(in) :: int_id, ispec, jspec
    character(path_len), allocatable, intent(out) :: fnames(:)
    integer :: interaction, nsshi, nsshj, nints, isorp, igrid, nzxi, nzxj
    character(2) :: auxisorp, auxzi, auxzj

    ! Retrieve basic info
    nints = twocenter_get_nints(int_id, ispec, jspec)
    nsshi = wf_atoms(ispec)%get_nshells()
    nsshj = wf_atoms(jspec)%get_nshells()
    nzxi = wf_atoms(ispec)%get_nz()
    nzxj = wf_atoms(jspec)%get_nz()
    write (auxzi,"(i2.2)") nzxi
    write (auxzj,"(i2.2)") nzxj
    interaction = ishft(1, int_id - 1)

    if (allocated(fnames)) deallocate (fnames)
    allocate (fnames(nints))
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_DENS_ATOM_SPH)
      do isorp = 1, nsshj
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp) = trim(TWOCENTER_ROOT(int_id))//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VNA_ATOM)
      do isorp = 0, nsshj
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(TWOCENTER_ROOT(int_id))//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_DENS_ONTOP, TWOCENTER_DENS_ONTOP_SPH)
      do isorp = 1, nsshi
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp) = trim(TWOCENTER_ROOT(int_id))//"l_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      do isorp = 1, nsshj
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + nsshi) = trim(TWOCENTER_ROOT(int_id))//"r_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VNA_ONTOP)
      do isorp = 0, nsshi
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(TWOCENTER_ROOT(int_id))//"l_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      do isorp = 0, nsshj
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + nsshi + 2) = trim(TWOCENTER_ROOT(int_id))//"r_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VXC)
      fnames(1) = trim(TWOCENTER_ROOT(int_id))//"_neutral."//auxzi//"."//auxzj//".dat"
      do isorp = 1, nsshi
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(TWOCENTER_ROOT(int_id))//"_ontopl_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      do isorp = 1, nsshj
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1 + nsshi) = trim(TWOCENTER_ROOT(int_id))//"_ontopr_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case default
      fnames(1) = trim(TWOCENTER_ROOT(int_id))//"."//auxzi//"."//auxzj//".dat"
    end select
  end subroutine twocenter_get_fnames

  subroutine twocenter_write_interaction(int_id, ispec, jspec, index_max, fnames, answer, names)
    integer, intent(in) :: int_id, ispec, jspec, index_max
    character(path_len), intent(in) :: fnames(:)
    real(dp), intent(in) :: answer(:,:,:)
    character(64), intent(in) :: names(:)
    integer :: interaction, issh, nsshi, nsshj, nints, isorp, igrid, nzxi, nzxj, index, io
    real(dp) :: rcuti, rcutj

    ! Retrieve basic info
    nints = twocenter_get_nints(int_id, ispec, jspec)
    nsshi = wf_atoms(ispec)%get_nshells()
    nsshj = wf_atoms(jspec)%get_nshells()
    rcuti = wf_atoms(ispec)%get_rcut()
    rcutj = wf_atoms(jspec)%get_rcut()
    nzxi = wf_atoms(ispec)%get_nz()
    nzxj = wf_atoms(jspec)%get_nz()
    interaction = ishft(1, int_id - 1)

    do isorp = 1, nints
      io = utils_open(fnames(isorp), "w")
      write (io, "(14x,i2,2x,i2,40x,'! Atomic numbers')") nzxi, nzxj
      write (io, "(2x,f14.6,2x,f14.6,28x,'! Cutoff radii')") rcuti, rcutj
      write (io, "(2x,f14.6,2x,f14.6,2x,i8,18x,'! zmin, zmax, nz')") 0.0_dp, rcuti + rcutj, TWOCENTER_NPOINTS_D
      if (interaction == TWOCENTER_VPP) then
        write (io, "(8x,i8,36x,'! npp')") pp_atoms(jspec)%get_nshells()
        write (io, "(1000ES16.8)") (pp_atoms(jspec)%get_cl(issh), issh = 1, nsshj)
      end if
      if (interaction == TWOCENTER_VXC) then
        write (io, "(8x,i8,2x,i8,26x,'! Number of shells')") nsshi, nsshj
        write (io, "(1000ES16.8)") (wf_atoms(ispec)%get_ref_charge(issh), issh = 1, nsshi)
        write (io, "(1000ES16.8)") (wf_atoms(jspec)%get_ref_charge(issh), issh = 1, nsshj)
      end if
      write (io, "('!')", advance="no")
      if (index_max == 0) then
        close (io)
        return
      end if
      do index = 1, index_max
        write (io, "(4x,a)", advance="no") trim(names(index))
      end do
      write (io, "(a)") ""
      do igrid = 1, TWOCENTER_NPOINTS_D
        write (io, "(1000ES16.8)") answer(isorp, :, igrid)
      end do
      close (io)
    end do
  end subroutine twocenter_write_interaction

  subroutine twocenter_calc(interactions, ispec, jspec)
    integer, intent(in) :: interactions, ispec, jspec
    integer :: interaction, int_id, isorp, nints, index_max
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    logical, allocatable :: exists(:)
    character(path_len), allocatable :: fnames(:)
    character(64), allocatable :: names(:)
    real(dp), allocatable :: answer(:,:,:)

    if (iand(TWOCENTER_NPOINTS_Z, 1) == 0 .or. iand(TWOCENTER_NPOINTS_RHO, 1) == 0) then
      write (stderr, "('[ERROR] ',(a))") "Number of integration points must be odd"
      error stop 1
    end if
    if (TWOCENTER_NPOINTS_Z == 1 .or. TWOCENTER_NPOINTS_RHO == 1) then
      write (stderr, "('[ERROR] ',(a))") "Number of integration points must be > 1"
      error stop 1
    end if

    call twocenter_get_interactions(interactions, twocenter_interactions)
    do int_id = 1, TWOCENTER_NUM_INTERACTIONS
      if (.not. twocenter_interactions(int_id)) cycle
      interaction = ishft(1, int_id - 1)
      call twocenter_get_fnames(int_id, ispec, jspec, fnames)
      nints = twocenter_get_nints(int_id, ispec, jspec)
      if (allocated(exists)) deallocate (exists)
      allocate (exists(nints))
      do isorp = 1, nints
        inquire (file=fnames(isorp), exist=exists(isorp))
      end do
      if (all(exists)) cycle
      do isorp = 1, nints
        write (stdout, "(2x,a)") "Computing "//trim(fnames(isorp))//"..."
      end do
      call twocenter_calc_interaction(int_id, ispec, jspec, index_max, answer, names)
      call twocenter_write_interaction(int_id, ispec, jspec, index_max, fnames, answer, names)
      write (stdout, "(1000a)", advance="no") repeat(achar(8)//achar(13), nints)
      do isorp = 1, nints
        write (stdout, "(2x,a)") "Computing "//trim(fnames(isorp))//"... Done!"
      end do
    end do ! int_id
    deallocate (exists, fnames, names, answer)
  end subroutine twocenter_calc

  real(dp) function cylindrical_harmonic(l, m, rho, z, r)
  ! s(sigma)  = 1
  ! p_sigma   = z/r
  ! p_pi      = rho/r
  ! d_sigma   = (2*z**2-rho**2)/r**2
  ! d_pi      = rho*z/r**2
  ! d_delta   = rho**2/r**2
  ! f_sigma   = z*(2*z**2-3*rho**2)/r**3
  ! f_pi      = rho*(4*z**2-rho**2)/r**3
  ! f_delta   = z*rho**2/r**3
  ! f_phi     = rho**3/r**3
    integer, intent(in) :: l, m
    real(dp), intent(in) :: rho, z, r

    if (r < tolerance) then
      cylindrical_harmonic = 0.0_dp
      return
    end if

    select case (l)
    case (0)
      cylindrical_harmonic = 1.0_dp
    case (1)
      select case (m)
      case (-1, 1)
        cylindrical_harmonic = 1.7320508075688772_dp*rho/r
      case (0)
        cylindrical_harmonic = 1.7320508075688772_dp*z/r
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (2)
      select case (m)
      case (-2, 2)
        cylindrical_harmonic = 1.9364916731037085_dp*(rho*rho)/(r*r)
      case (-1, 1)
        cylindrical_harmonic = 3.872983346207417_dp*(rho*z)/(r*r)
      case (0)
        cylindrical_harmonic = 2.23606797749979_dp*(2.0_dp*z*z - rho*rho)/(r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (3)
      select case (m)
      case (-3, 3)
        cylindrical_harmonic = 2.091650066335189_dp*(rho*rho*rho)/(r*r*r)
      case (-2, 2)
        cylindrical_harmonic = 5.123475382979799_dp*(rho*rho*z)/(r*r*r)
      case (-1, 1)
        cylindrical_harmonic = 0.6123724356957945_dp*rho*(4.0_dp*z*z - rho*rho)/(r*r*r)
      case (0)
        cylindrical_harmonic = 1.3228756555322954_dp*z*(2.0_dp*z*z - 3.0_dp*rho*rho*rho)/(r*r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
      case default
        cylindrical_harmonic = 0.0_dp
    end select
    cylindrical_harmonic = cylindrical_harmonic * sqinv4pi
  end function cylindrical_harmonic
end module twocenter
