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
    &                   INDICES_TWOCENTER_DIPY, INDICES_TWOCENTER_COULOMB, INDICES_TWOCENTER_SPH
  use :: xc, only: xc_calc, xc_isgga
  use :: wavefunctions, only: wf_atoms, wf_dens
  use :: potentials, only: pot_atoms
  use :: pseudopotentials, only: pp_atoms
  implicit none
  private
  public :: twocenter_calc

  integer, parameter, public :: TWOCENTER_NPOINTS_D = 107
  integer, parameter, public :: TWOCENTER_NPOINTS_Z = 213
  integer, parameter, public :: TWOCENTER_NPOINTS_RHO = 107

  integer, parameter :: TWOCENTER_NUM_INTERACTIONS = 11
  integer, parameter, public :: TWOCENTER_DENS_ATOM  = ishft(1, 0)
  integer, parameter, public :: TWOCENTER_DENS_ONTOP = ishft(1, 1)
  integer, parameter, public :: TWOCENTER_OVERLAP    = ishft(1, 2)
  integer, parameter, public :: TWOCENTER_VNA_ATOM   = ishft(1, 3)
  integer, parameter, public :: TWOCENTER_VNA_ONTOP  = ishft(1, 4)
  integer, parameter, public :: TWOCENTER_VPP        = ishft(1, 5)
  integer, parameter, public :: TWOCENTER_VXC        = ishft(1, 6)
  integer, parameter, public :: TWOCENTER_DIP_Z      = ishft(1, 7)
  integer, parameter, public :: TWOCENTER_DIP_X      = ishft(1, 8)
  integer, parameter, public :: TWOCENTER_DIP_Y      = ishft(1, 9)
  integer, parameter, public :: TWOCENTER_COULOMB    = ishft(1, 10)
  character(32), parameter :: TWOCENTER_ROOT(TWOCENTER_NUM_INTERACTIONS) = [ &
    &  "den_atom                        ", &
    &  "den_ontop                       ", &
    &  "overlap                         ", &
    &  "vna_atom                        ", &
    &  "vna_ontop                       ", &
    &  "vnl                             ", &
    &  "vxc                             ", &
    &  "dipole_z                        ", &
    &  "dipole_x                        ", &
    &  "dipole_y                        ", &
    &  "coulomb                         " &
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
    integer :: interaction, nssh1, nssh2
    interaction = ishft(1, int_id - 1)
    nssh1 = wf_atoms(ispec)%get_nshells()
    nssh2 = wf_atoms(jspec)%get_nshells()
    select case (interaction)
    case (TWOCENTER_DENS_ATOM)
      twocenter_get_nints = nssh2
    case (TWOCENTER_VNA_ATOM)
      twocenter_get_nints = nssh2 + 1
    case (TWOCENTER_DENS_ONTOP)
      twocenter_get_nints = nssh1 + nssh2
    case (TWOCENTER_VXC)
      twocenter_get_nints = nssh1 + nssh2 + 1
    case (TWOCENTER_VNA_ONTOP)
      twocenter_get_nints = nssh1 + nssh2 + 2
    case default
      twocenter_get_nints = 1
    end select
  end function twocenter_get_nints

  pure integer function twocenter_get_index_type(int_id, issph)
    integer, intent(in) :: int_id
    logical, intent(in) :: issph
    integer :: interaction
    if (issph) then
      twocenter_get_index_type = INDICES_TWOCENTER_SPH
      return
    end if
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

  subroutine twocenter_calc_interaction(int_id, ispec, jspec, issph, index_max, answer, names)
    integer, intent(in) :: int_id, ispec, jspec
    logical, intent(in) :: issph
    integer, intent(out) :: index_max
    real(dp), allocatable, intent(out) :: answer(:,:,:)
    character(64), allocatable, intent(out) :: names(:)
    integer :: nssh1, nssh2, nints, interaction, igrid, iz, irho, isorp, idx, nz, nrho
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    real(dp) :: d, dd, dmax, rcut1, rcut2, zmin, zmax, dz, drho, z1, z2, z12, z22, r1, r2, rhomult, rhomax, &
      &         zmult, factor, psi1, psi2, cyl1, cyl2, psimult, tmp, rho, rho2, &
      &         dens, exc, vxc, dexcrho, dvxcrho, dexcsigma, dvxcsigma
    real(dp) :: ddens(2), dddens(2, 2)
    integer, allocatable :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:), ls1(:), ls2(:)
    real(dp), allocatable :: fofr(:)

    interaction = ishft(1, int_id - 1)

    ! Retrieve basic info
    nssh1 = wf_atoms(ispec)%get_nshells()
    nssh2 = wf_atoms(jspec)%get_nshells()
    rcut1 = wf_atoms(ispec)%get_rcut()
    ls1 = [(wf_atoms(ispec)%get_angular_momentum(idx), idx = 1, wf_atoms(ispec)%get_nshells())]
    select case (interaction)
    case (TWOCENTER_VPP)
      rcut2 = pp_atoms(jspec)%get_rcut()
    case default
      rcut2 = wf_atoms(jspec)%get_rcut()
    end select
    select case (interaction)
    case (TWOCENTER_DENS_ATOM, TWOCENTER_VNA_ATOM)
      ls2 = [(wf_atoms(ispec)%get_angular_momentum(idx), idx = 1, wf_atoms(ispec)%get_nshells())]
    case (TWOCENTER_VPP)
      ls2 = [(pp_atoms(jspec)%get_angular_momentum(idx), idx = 1, pp_atoms(jspec)%get_nshells())]
    case default
      ls2 = [(wf_atoms(jspec)%get_angular_momentum(idx), idx = 1, wf_atoms(jspec)%get_nshells())]
    end select

    ! Set dimensions
    nints = twocenter_get_nints(int_id, ispec, jspec)
    call indices_twocenter_set(twocenter_get_index_type(int_id, issph), ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names=names)
    deallocate (ls1, ls2)
    if (index_max == 0) return
    if (allocated(answer)) deallocate (answer)
    allocate (fofr(nints), answer(nints, index_max, TWOCENTER_NPOINTS_D))
    answer = 0.0_dp

    ! Prepare integral
    dmax = rcut1 + rcut2
    dd = dmax/real(TWOCENTER_NPOINTS_D - 1, kind=dp)
    drho = min(rcut1, rcut2)/real(TWOCENTER_NPOINTS_RHO - 1, kind=dp)

    do igrid = 1, TWOCENTER_NPOINTS_D
      d = real(igrid - 1, kind=dp)*dd

      select case (interaction)
      case (TWOCENTER_COULOMB)
        zmin = -rcut1
        zmax = rcut1
        rhomax = min(rcut1, rcut2)
      case default
        zmin = max(-rcut1, d - rcut2)
        zmax = min(rcut1, d + rcut2)
        rhomax = rcut1
      end select

      ! Strictly define what the density of the mesh should be. Make the density of
      ! the number of points equivalent for all cases. Change the number of points
      ! to be integrated to be dependent upon the distance between the centers and
      ! this defined density.
      dz = (rcut1 + rcut2)/real(TWOCENTER_NPOINTS_Z - 1, kind=dp)
      drho = max(rcut1, rcut2)/real(TWOCENTER_NPOINTS_RHO - 1, kind=dp)
      nz = int((zmax - zmin)/dz) + 1
      nrho = int(rhomax/drho) + 1
      if (iand(nz, 1) == 0) then
        nz = nz + 1
        dz = (zmax - zmin)/real(nz - 1, kind=dp)
      end if
      if (iand(nrho, 1) == 0) then
        nrho = nrho + 1
        drho = rhomax/real(nrho - 1, kind=dp)
      end if

      do iz = 1, nz
        z1 = zmin + real(iz - 1, kind=dp)*dz
        z2 = z1 - d
        z12 = z1*z1
        z22 = z2*z2
        zmult = 0.66666666666666666667_dp*dz
        if (iz == 1 .or. iz == nz) then
          zmult = 0.5_dp*zmult
        else if (iand(iz, 1) == 0) then
          zmult = 2.0_dp*zmult
        end if

        do irho = 2, nrho ! We can ignore 0 because of factor
          rho = real(irho - 1, kind=dp)*drho
          rhomult = 0.66666666666666666667_dp*drho
          if (irho == nrho) then
            rhomult = 0.5_dp*rhomult
          else if (iand(irho, 1) == 0) then
            rhomult = 2.0_dp*rhomult
          end if
          factor = zmult*rhomult*twopi

          rho2 = rho*rho
          r1 = sqrt(z12 + rho2)
          if (r1 >= rcut1) cycle
          r2 = sqrt(z22 + rho2)

          ! The integral in all its glory
          select case (interaction)
          case (TWOCENTER_DENS_ATOM)
            do isorp = 1, nssh2
              tmp = wf_atoms(jspec)%get_psi(isorp, r2)
              fofr(isorp) = inv4pi*tmp*tmp
            end do
          case (TWOCENTER_DENS_ONTOP)
            do isorp = 1, nssh1
              tmp = wf_atoms(ispec)%get_psi(isorp, r1)
              fofr(isorp) = inv4pi*tmp*tmp
            end do
            do isorp = 1, nssh2
              tmp = wf_atoms(jspec)%get_psi(isorp, r2)
              fofr(isorp + nssh1) = inv4pi*tmp*tmp
            end do
          case (TWOCENTER_VNA_ATOM)
            do isorp = 0, nssh2
              fofr(isorp + 1) = pot_atoms(jspec)%get_vnn(isorp, r2)
            end do
          case (TWOCENTER_VNA_ONTOP)
            do isorp = 0, nssh1
              fofr(isorp + 1) = pot_atoms(ispec)%get_vnn(isorp, r1)
            end do
            do isorp = 0, nssh2
              fofr(isorp + nssh1 + 2) = pot_atoms(jspec)%get_vnn(isorp, r2)
            end do
          case (TWOCENTER_VXC)
            if (xc_isgga()) then
              call wf_dens(ispec, jspec, rho, z1, z2, r1, r2, dens, ddens, dddens)
              call xc_calc(dens, ddens, dddens, exc, vxc, dexcrho, dexcsigma, dvxcrho, dvxcsigma)
            else
              call wf_dens(ispec, jspec, rho, z1, z2, r1, r2, dens)
              call xc_calc(dens, exc, vxc, dexcrho, dvxcrho)
              dexcsigma = 0.0_dp
              dvxcsigma = 0.0_dp
            end if
            fofr(1) = vxc
            do isorp = 1, nssh1
              tmp = wf_atoms(ispec)%get_psi(isorp, r1)
              fofr(isorp + 1) = inv4pi*tmp*tmp*dvxcrho
              if (r1 > tolerance) fofr(isorp + 1) = fofr(isorp + 1) + dvxcsigma*4.0_dp*inv4pi * &
                  &               tmp*wf_atoms(ispec)%get_psi(isorp, r1, order=1)*(ddens(1)*rho + ddens(2)*z1)/r1
            end do
            do isorp = 1, nssh2
              tmp = wf_atoms(jspec)%get_psi(isorp, r2)
              fofr(isorp + 1 + nssh1) = inv4pi*tmp*tmp*dvxcrho
              if (r2 > tolerance) fofr(isorp + 1 + nssh1) = fofr(isorp + 1 + nssh1) + dvxcsigma*4.0_dp*inv4pi * &
                  &               tmp*wf_atoms(jspec)%get_psi(isorp, r2, order=1)*(ddens(1)*rho + ddens(2)*z2)/r2
            end do
          case (TWOCENTER_DIP_Z)
            fofr(1) = z1 - 0.5_dp*d
          case (TWOCENTER_DIP_X, TWOCENTER_DIP_Y)
            fofr(1) = rho
          case default
            fofr(1) = 1.0_dp
          end select

          do idx = 1, index_max
            psi1 = wf_atoms(ispec)%get_psi(s1(idx), r1)
            select case (interaction)
            case (TWOCENTER_DENS_ATOM, TWOCENTER_VNA_ATOM)
              psi2 = wf_atoms(ispec)%get_psi(s2(idx), r1)
            case (TWOCENTER_VPP)
              psi2 = pp_atoms(jspec)%get_vpp(s2(idx), r2)
            case (TWOCENTER_COULOMB)
              psi2 = psi1*pot_atoms(jspec)%get_vnn(s2(idx), r2)
            case default
              psi2 = wf_atoms(jspec)%get_psi(s2(idx), r2)
            end select

            cyl1 = cylindrical_harmonic(int_id, l1(idx), m1(idx), rho, z1, r1)
            cyl2 = cylindrical_harmonic(int_id, l2(idx), m2(idx), rho, z2, r2)
            psimult = psi1*psi2
            if (issph .and. (psimult < 0)) psimult = -psimult
            psimult = psimult*cyl1*cyl2
            answer(:, idx, igrid) = answer(:, idx, igrid) + fofr(:)*psimult*factor*rho
          end do ! idx
        end do ! irho
      end do ! iz
    end do ! igrid
    deallocate (s1, s2, l1, l2, m1, m2, fofr)
  end subroutine twocenter_calc_interaction

  subroutine twocenter_get_fnames(int_id, ispec, jspec, issph, fnames)
    integer, intent(in) :: int_id, ispec, jspec
    logical, intent(in) :: issph
    character(path_len), allocatable, intent(out) :: fnames(:)
    integer :: interaction, nssh1, nssh2, nints, isorp, igrid, nzxi, nzxj
    character(1) :: csph
    character(2) :: auxisorp, auxzi, auxzj
    character(32) :: root

    ! Retrieve basic info
    nints = twocenter_get_nints(int_id, ispec, jspec)
    nssh1 = wf_atoms(ispec)%get_nshells()
    nssh2 = wf_atoms(jspec)%get_nshells()
    nzxi = wf_atoms(ispec)%get_nz()
    nzxj = wf_atoms(jspec)%get_nz()
    write (auxzi,"(i2.2)") nzxi
    write (auxzj,"(i2.2)") nzxj
    interaction = ishft(1, int_id - 1)

    csph = " "
    if (issph) csph = "S"
    if (allocated(fnames)) deallocate (fnames)
    allocate (fnames(nints))
    select case (interaction)
    case (TWOCENTER_DENS_ATOM)
      root = trim(TWOCENTER_ROOT(int_id))//trim(csph)
      do isorp = 1, nssh2
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VNA_ATOM)
      root = trim(TWOCENTER_ROOT(int_id))//trim(csph)
      do isorp = 0, nssh2
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_DENS_ONTOP)
      root = trim(TWOCENTER_ROOT(int_id))//"l"//trim(csph)
      do isorp = 1, nssh1
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      root = trim(TWOCENTER_ROOT(int_id))//"r"//trim(csph)
      do isorp = 1, nssh2
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + nssh1) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VNA_ONTOP)
      root = trim(TWOCENTER_ROOT(int_id))//"l"//trim(csph)
      do isorp = 0, nssh1
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      root = trim(TWOCENTER_ROOT(int_id))//"r"//trim(csph)
      do isorp = 0, nssh2
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + nssh1 + 2) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case (TWOCENTER_VXC)
      root = trim(TWOCENTER_ROOT(int_id))//"_neutral"//trim(csph)
      fnames(1) = trim(root)//"."//auxzi//"."//auxzj//".dat"
      root = trim(TWOCENTER_ROOT(int_id))//"_ontopl"//trim(csph)
      do isorp = 1, nssh1
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
      root = trim(TWOCENTER_ROOT(int_id))//"_ontopr"//trim(csph)
      do isorp = 1, nssh2
        write (auxisorp,"(i2.2)") isorp
        fnames(isorp + 1 + nssh1) = trim(root)//"_"//auxisorp//"."//auxzi//"."//auxzj//".dat"
      end do
    case default
      root = trim(TWOCENTER_ROOT(int_id))//trim(csph)
      fnames(1) = trim(root)//"."//auxzi//"."//auxzj//".dat"
    end select
  end subroutine twocenter_get_fnames

  subroutine twocenter_write_interaction(int_id, ispec, jspec, index_max, fnames, answer, names)
    integer, intent(in) :: int_id, ispec, jspec, index_max
    character(path_len), intent(in) :: fnames(:)
    real(dp), intent(in) :: answer(:,:,:)
    character(64), intent(in) :: names(:)
    integer :: interaction, issh, nssh1, nssh2, nints, isorp, igrid, nzxi, nzxj, index, io
    real(dp) :: rcut1, rcut2

    ! Retrieve basic info
    nints = twocenter_get_nints(int_id, ispec, jspec)
    nssh1 = wf_atoms(ispec)%get_nshells()
    nssh2 = wf_atoms(jspec)%get_nshells()
    rcut1 = wf_atoms(ispec)%get_rcut()
    rcut2 = wf_atoms(jspec)%get_rcut()
    nzxi = wf_atoms(ispec)%get_nz()
    nzxj = wf_atoms(jspec)%get_nz()
    interaction = ishft(1, int_id - 1)
    if (interaction == TWOCENTER_VPP) then
      nssh2 = pp_atoms(jspec)%get_nshells()
      rcut2 = pp_atoms(jspec)%get_rcut()
    end if

    do isorp = 1, nints
      io = utils_open(fnames(isorp), "w")
      write (io, "(14x,i2,2x,i2,40x,'! Atomic numbers')") nzxi, nzxj
      write (io, "(2x,f14.6,2x,f14.6,28x,'! Cutoff radii')") rcut1, rcut2
      write (io, "(2x,f14.6,2x,f14.6,2x,i8,18x,'! zmin, zmax, nz')") 0.0_dp, rcut1 + rcut2, TWOCENTER_NPOINTS_D
      if (interaction == TWOCENTER_VPP) then
        write (io, "(8x,i8,36x,'! npp')") nssh2
        write (io, "(1000ES16.8)") (pp_atoms(jspec)%get_cl(issh), issh = 1, nssh2)
      end if
      if (interaction == TWOCENTER_VXC) then
        write (io, "(8x,i8,2x,i8,26x,'! Number of shells')") nssh1, nssh2
        write (io, "(1000ES16.8)", advance="no") (wf_atoms(ispec)%get_ref_charge(issh), issh = 1, nssh1)
        write (io, "(a)") " ! Reference charges atom 1"
        write (io, "(1000ES16.8)", advance="no") (wf_atoms(jspec)%get_ref_charge(issh), issh = 1, nssh2)
        write (io, "(a)") " ! Reference charges atom 2"
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

  subroutine twocenter_calc(interactions, ispec, jspec, is_spherical)
    integer, intent(in) :: interactions, ispec, jspec
    logical, intent(in), optional :: is_spherical
    integer :: interaction, int_id, isorp, nints, index_max
    logical :: issph
    logical :: twocenter_interactions(TWOCENTER_NUM_INTERACTIONS)
    logical, allocatable :: exists(:)
    character(path_len), allocatable :: fnames(:)
    character(64), allocatable :: names(:)
    real(dp), allocatable :: answer(:,:,:)

    issph = .false.
    if (present(is_spherical)) issph = is_spherical

    call twocenter_get_interactions(interactions, twocenter_interactions)
    do int_id = 1, TWOCENTER_NUM_INTERACTIONS
      if (.not. twocenter_interactions(int_id)) cycle
      interaction = ishft(1, int_id - 1)
      call twocenter_get_fnames(int_id, ispec, jspec, issph, fnames)
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
      call twocenter_calc_interaction(int_id, ispec, jspec, issph, index_max, answer, names)
      call twocenter_write_interaction(int_id, ispec, jspec, index_max, fnames, answer, names)
      write (stdout, "(1000a)", advance="no") repeat(achar(8)//achar(13), nints)
      do isorp = 1, nints
        write (stdout, "(2x,a)") "Computing "//trim(fnames(isorp))//"... Done!"
      end do
    end do ! int_id
    deallocate (exists, fnames, names, answer)
  end subroutine twocenter_calc

  real(dp) function cylindrical_harmonic(int_id, l, m, rho, z, r)
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
    integer, intent(in) :: int_id, l, m
    real(dp), intent(in) :: rho, z, r
    integer :: interaction
    interaction = ishft(1, int_id - 1)

    if (interaction == TWOCENTER_COULOMB) then
      cylindrical_harmonic = sqinv4pi
      return
    end if

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
        cylindrical_harmonic = 1.224744871391589_dp*rho/r
      case (0)
        cylindrical_harmonic = 1.7320508075688772_dp*z/r
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (2)
      select case (m)
      case (-2, 2)
        cylindrical_harmonic = 1.3693063937629153_dp*(rho*rho)/(r*r)
      case (-1, 1)
        cylindrical_harmonic = 2.7386127875258306_dp*(rho*z)/(r*r)
      case (0)
        cylindrical_harmonic = 1.118033988749895_dp*(2.0_dp*z*z - rho*rho)/(r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
    case (3)
      select case (m)
      case (-3, 3)
        cylindrical_harmonic = 1.479019945774904_dp*(rho*rho*rho)/(r*r*r)
      case (-2, 2)
        cylindrical_harmonic = 3.6228441865473595_dp*(rho*rho*z)/(r*r*r)
      case (-1, 1)
        cylindrical_harmonic = 1.14564392373896_dp*rho*(4.0_dp*z*z - rho*rho)/(r*r*r)
      case (0)
        cylindrical_harmonic = 1.3228756555322954_dp*z*(2.0_dp*z*z - 3.0_dp*rho*rho*rho)/(r*r*r)
      case default
        cylindrical_harmonic = 0.0_dp
      end select
      case default
        cylindrical_harmonic = 0.0_dp
    end select
    cylindrical_harmonic = cylindrical_harmonic * sqinv4pi

    ! Some signs and factors need to be adjusted for dipX and dipY
    if (interaction == TWOCENTER_DIP_X .or. interaction == TWOCENTER_DIP_Y) then
      if (iand(m, 1) == 0) cylindrical_harmonic = cylindrical_harmonic*invsq2
      if (interaction == TWOCENTER_DIP_Y) then
        if (l == 2 .and. m == 2 .or. &
          & l == 3 .and. (m == -3 .or. m == 2 .or. m == 3)) cylindrical_harmonic = -cylindrical_harmonic
      end if
    end if
  end function cylindrical_harmonic
end module twocenter
