! copyright info:
!
!                             @Copyright 2026
!                           fireball-qmd project
!
! fireball-qmd is a free (GPLv3) open project.
!
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

! input.f90
! Module Description
! ===========================================================================
!       Reads the single input file that 'fireballpy.basis' writes for
! 'create', and hands back everything create needs in one structure.
!
!       This replaces the old pair of files (an index 'create.input' listing
! one '<El>.input' per species, both positional) and the three independent
! parsers that used to walk them -- wf_init, pot_init and pp_init each
! reopened 'create.input' and skipped a different number of lines, so adding
! a field meant editing three hardcoded counters at once, and forgetting one
! did not fail: it silently read the neighbouring field.
!
!       The format is sections plus 'key = value', so order inside a section
! does not matter and unknown keys are reported instead of shifting
! everything by one line.
!
! ===========================================================================
!
! Module Declaration
! ===========================================================================
module input
  use, intrinsic :: iso_fortran_env, only: dp => real64, iostat_end
  use :: utils, only: utils_open, utils_itos, utils_log, UTILS_LOG_ERROR
  implicit none
  private

  public :: input_t, input_species_t, input_shell_t, input_read

  integer, parameter :: INPUT_MAXLINE = 1024
  integer, parameter :: INPUT_MAXTOK = 64
  integer, parameter :: INPUT_TOKLEN = 256
  integer, parameter :: INPUT_NAMELEN = 32
  ! Bump this when the format changes in a way old readers cannot handle.
  integer, parameter, public :: INPUT_VERSION = 1

  type :: input_shell_t
    integer :: l = -1
    ! Reference ("neutral") charge of the shell. It is NOT the occupation the
    ! orbital was generated with: that one is a 'begin' artefact used to
    ! compress the basis and has no meaning here.
    real(dp) :: qref = 0.0_dp
    character(:), allocatable :: fwf, fna
  end type input_shell_t

  type :: input_species_t
    character(2) :: symbol = '  '
    integer :: znuc = 0
    integer :: nshells = 0
    character(:), allocatable :: ppfile, na0file
    type(input_shell_t), allocatable :: shells(:)
  end type input_species_t

  type :: input_t
    integer :: version = 0
    character(:), allocatable :: xc_names(:)
    real(dp), allocatable :: xc_weights(:)
    type(input_species_t), allocatable :: species(:)
  end type input_t

contains

  ! =========================================================================
  ! Read the whole file and return it parsed. Any problem is fatal: a basis
  ! built from a half-understood input is worse than no basis.
  ! =========================================================================
  type(input_t) function input_read(path) result(inp)
    character(*), intent(in) :: path
    integer :: io, ios, i, isp, ish, ntok
    character(INPUT_MAXLINE) :: line
    character(INPUT_TOKLEN) :: section, key, value
    character(INPUT_TOKLEN) :: toks(INPUT_MAXTOK)
    logical :: has_count

    io = utils_open(path, "r")
    section = ''
    isp = 0
    ish = 0
    has_count = .false.

    do
      read (io, '(a)', iostat=ios) line
      if (ios == iostat_end) exit
      if (ios /= 0) call fail("could not read "//trim(path))
      call strip(line)
      if (len_trim(line) == 0) cycle

      ! --- section header -------------------------------------------------
      if (line(1:1) == '[') then
        i = index(line, ']')
        if (i < 3) call fail("malformed section: "//trim(line))
        section = line(2:i - 1)
        if (section(1:8) == 'species.') then
          if (.not. has_count) call fail("[species] count must come before [species.N]")
          read (section(9:), *, iostat=ios) isp
          if (ios /= 0) call fail("malformed species index: "//trim(section))
          if (isp < 1 .or. isp > size(inp%species)) call fail("species index out of range: "//trim(section))
          ish = 0
        end if
        cycle
      end if

      ! --- key = value ----------------------------------------------------
      i = index(line, '=')
      if (i == 0) call fail("expected 'key = value', got: "//trim(line))
      key = adjustl(line(1:i - 1))
      value = adjustl(line(i + 1:))

      select case (trim(section))

      case ('meta')
        select case (trim(key))
        case ('version')
          read (value, *, iostat=ios) inp%version
          if (ios /= 0) call fail("bad version: "//trim(value))
          if (inp%version /= INPUT_VERSION) call fail("unsupported input version "//trim(value))
        case ('generator')
          ! informative only
        case default
          call fail("unknown key in [meta]: "//trim(key))
        end select

      case ('xc')
        call tokenize(value, toks, ntok)
        if (ntok == 0) call fail("empty value for "//trim(key)//" in [xc]")
        select case (trim(key))
        case ('names')
          allocate (character(INPUT_NAMELEN) :: inp%xc_names(ntok))
          do i = 1, ntok
            inp%xc_names(i) = trim(toks(i))
          end do
        case ('weights')
          allocate (inp%xc_weights(ntok))
          do i = 1, ntok
            read (toks(i), *, iostat=ios) inp%xc_weights(i)
            if (ios /= 0) call fail("bad weight: "//trim(toks(i)))
          end do
        case default
          call fail("unknown key in [xc]: "//trim(key))
        end select

      case ('species')
        if (trim(key) /= 'count') call fail("unknown key in [species]: "//trim(key))
        read (value, *, iostat=ios) i
        if (ios /= 0 .or. i < 1) call fail("bad species count: "//trim(value))
        allocate (inp%species(i))
        has_count = .true.

      case default
        if (section(1:8) /= 'species.') call fail("unknown section: ["//trim(section)//"]")
        call read_species_key(inp%species(isp), key, value, ish)

      end select
    end do
    close (io)

    call check(inp, path)
  end function input_read

  ! =========================================================================
  subroutine read_species_key(sp, key, value, ish)
    type(input_species_t), intent(inout) :: sp
    character(*), intent(in) :: key, value
    integer, intent(inout) :: ish
    integer :: ios, ntok
    character(INPUT_TOKLEN) :: toks(INPUT_MAXTOK)

    select case (trim(key))
    case ('symbol')
      sp%symbol = trim(value)
    case ('znuc')
      read (value, *, iostat=ios) sp%znuc
      if (ios /= 0) call fail("bad znuc: "//trim(value))
    case ('ppfile')
      sp%ppfile = trim(value)
    case ('na0file')
      sp%na0file = trim(value)
    case ('nshells')
      read (value, *, iostat=ios) sp%nshells
      if (ios /= 0 .or. sp%nshells < 1) call fail("bad nshells: "//trim(value))
      allocate (sp%shells(sp%nshells))
    case ('shell')
      ! shell = <l> <qref> <wavefunction file> <neutral potential file>
      if (.not. allocated(sp%shells)) call fail("'nshells' must come before 'shell'")
      ish = ish + 1
      if (ish > sp%nshells) call fail("more 'shell' lines than 'nshells'")
      call tokenize(value, toks, ntok)
      if (ntok /= 4) call fail("'shell' needs 4 fields (l qref wf na), got: "//trim(value))
      read (toks(1), *, iostat=ios) sp%shells(ish)%l
      if (ios /= 0) call fail("bad l: "//trim(toks(1)))
      read (toks(2), *, iostat=ios) sp%shells(ish)%qref
      if (ios /= 0) call fail("bad reference charge: "//trim(toks(2)))
      sp%shells(ish)%fwf = trim(toks(3))
      sp%shells(ish)%fna = trim(toks(4))
    case default
      call fail("unknown key in [species.N]: "//trim(key))
    end select
  end subroutine read_species_key

  ! =========================================================================
  ! Everything that must be there, in one place, so a missing field is a
  ! clear message and not a crash three modules later.
  ! =========================================================================
  subroutine check(inp, path)
    type(input_t), intent(in) :: inp
    character(*), intent(in) :: path
    integer :: i, j

    if (inp%version /= INPUT_VERSION) call fail(trim(path)//": missing [meta] version")
    if (.not. allocated(inp%xc_names)) call fail(trim(path)//": missing [xc] names")
    if (.not. allocated(inp%xc_weights)) call fail(trim(path)//": missing [xc] weights")
    if (size(inp%xc_names) /= size(inp%xc_weights)) &
      call fail(trim(path)//": [xc] names and weights differ in length")
    if (.not. allocated(inp%species)) call fail(trim(path)//": missing [species] count")

    do i = 1, size(inp%species)
      associate (sp => inp%species(i))
        if (sp%znuc < 1) call fail(trim(path)//": species "//utils_itos(i)//" has no znuc")
        if (sp%symbol == '  ') call fail(trim(path)//": species "//utils_itos(i)//" has no symbol")
        if (.not. allocated(sp%ppfile)) call fail(trim(path)//": species "//utils_itos(i)//" has no ppfile")
        if (.not. allocated(sp%na0file)) call fail(trim(path)//": species "//utils_itos(i)//" has no na0file")
        if (.not. allocated(sp%shells)) call fail(trim(path)//": species "//utils_itos(i)//" has no shells")
        do j = 1, sp%nshells
          if (sp%shells(j)%l < 0) call fail(trim(path)//": species "//utils_itos(i)//" is missing shell "//utils_itos(j))
        end do
      end associate
    end do
  end subroutine check

  ! =========================================================================
  ! Drop a trailing '#' comment and any tabs, so the file can be annotated.
  ! =========================================================================
  subroutine strip(line)
    character(*), intent(inout) :: line
    integer :: i
    do i = 1, len(line)
      if (iachar(line(i:i)) == 9) line(i:i) = ' '
    end do
    i = index(line, '#')
    if (i > 0) line(i:) = ' '
    line = adjustl(line)
  end subroutine strip

  ! =========================================================================
  subroutine tokenize(str, toks, ntok)
    character(*), intent(in) :: str
    character(INPUT_TOKLEN), intent(out) :: toks(:)
    integer, intent(out) :: ntok
    integer :: i, ini
    logical :: inside

    ntok = 0
    inside = .false.
    ini = 1
    do i = 1, len_trim(str) + 1
      if (i > len_trim(str) .or. str(i:i) == ' ') then
        if (inside) then
          ntok = ntok + 1
          if (ntok > size(toks)) call fail("too many fields in: "//trim(str))
          toks(ntok) = str(ini:i - 1)
          inside = .false.
        end if
      else
        if (.not. inside) then
          ini = i
          inside = .true.
        end if
      end if
    end do
  end subroutine tokenize

  ! =========================================================================
  subroutine fail(msg)
    character(*), intent(in) :: msg
    call utils_log(UTILS_LOG_ERROR, msg)
    error stop 1
  end subroutine fail

end module input
