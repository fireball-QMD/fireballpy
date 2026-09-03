module indices
  use :: constants, only: snames, pnames, dnames, fnames
  implicit none
  private
  public :: indices_onecenter_set, indices_twocenter_set

  integer, parameter, public :: INDICES_TWOCENTER = 0
  integer, parameter, public :: INDICES_TWOCENTER_DIPX = 1
  integer, parameter, public :: INDICES_TWOCENTER_DIPY = 2
  integer, parameter, public :: INDICES_TWOCENTER_COULOMB = 3
  integer, parameter, public :: INDICES_TWOCENTER_SPH = 4

contains

  pure subroutine indices_onecenter_set(ls, index_max, s1, s2, l12, names)
    integer, intent(in) :: ls(:)
    integer, intent(out) :: index_max
    integer, allocatable, intent(out) :: s1(:), s2(:), l12(:)
    character(64), allocatable, intent(out), optional :: names(:)
    integer :: nssh, issh, jssh, ix, l1tmp, l2tmp
    integer :: buffer(3, 1024)
    character(64) :: auxname

    if (allocated(s1)) deallocate (s1)
    if (allocated(s2)) deallocate (s2)
    if (allocated(l12)) deallocate (l12)
    if (present(names)) then
      if (allocated(names)) deallocate (names)
    end if
    nssh = size(ls)

    index_max = 0
    do issh = 1, nssh
      l1tmp = ls(issh)
      do jssh = 1, nssh
        l2tmp = ls(jssh)
        if (l1tmp /= l2tmp) cycle
        index_max = index_max + 1
        buffer(1, index_max) = issh
        buffer(2, index_max) = jssh
        buffer(3, index_max) = l1tmp
      end do
    end do
    if (index_max == 0) return
    allocate (s1(index_max), s2(index_max), l12(index_max))
    do ix = 1, index_max
      s1(ix) = buffer(1, ix)
      s2(ix) = buffer(2, ix)
      l12(ix) = buffer(3, ix)
    end do
    if (present(names)) then
      allocate (names(index_max))
      do ix = 1, index_max
        auxname = indices_get_name(l12(ix), 0, omit_suffix=.true.)
        names(ix) = trim(auxname)
      end do
    end if
  end subroutine indices_onecenter_set

  pure subroutine indices_twocenter_set(index_type, ls1, ls2, index_max, s1, s2, l1, l2, m1, m2, names)
    integer, intent(in) :: index_type
    integer, intent(in) :: ls1(:), ls2(:)
    integer, intent(out) :: index_max
    integer, allocatable, intent(out) :: s1(:), s2(:), l1(:), l2(:), m1(:), m2(:)
    character(64), allocatable, intent(out), optional :: names(:)
    integer :: nssh1, nssh2, issh, jssh, ix, l1tmp, l2tmp, m1tmp, m2tmp
    integer :: buffer(6, 1024)
    character(64) :: auxname1, auxname2

    if (allocated(s1)) deallocate (s1)
    if (allocated(s2)) deallocate (s2)
    if (allocated(l1)) deallocate (l1)
    if (allocated(l2)) deallocate (l2)
    if (allocated(m1)) deallocate (m1)
    if (allocated(m2)) deallocate (m2)
    if (present(names)) then
      if (allocated(names)) deallocate (names)
    end if
    nssh1 = size(ls1)
    nssh2 = size(ls2)

    ! Spherical is special
    if (index_type == INDICES_TWOCENTER_SPH) then
      index_max = nssh1*nssh2
      if (index_max == 0) return
      allocate (s1(index_max), s2(index_max), l1(index_max), l2(index_max), m1(index_max), m2(index_max))
      if (present(names)) allocate (names(index_max))
      ix = 0
      do issh = 1, nssh1
        do jssh = 1, nssh2
          ix = ix + 1
          if (present(names)) then
            auxname1 = indices_get_name(ls1(issh), 0, omit_suffix=.true.)
            auxname2 = indices_get_name(ls2(jssh), 0, omit_suffix=.true.)
            names(ix) = trim(auxname1)//trim(auxname2)
          end if
          s1(ix) = issh
          l1(ix) = 0
          m1(ix) = 0
          s2(ix) = jssh
          l2(ix) = 0
          m2(ix) = 0
        end do
      end do
      return
    end if

    index_max = 0
    do issh = 1, nssh1
      l1tmp = ls1(issh)
      do jssh = 1, nssh2
        l2tmp = ls2(jssh)
        do m1tmp = -l1tmp, l1tmp
          do m2tmp = -l2tmp, l2tmp
            select case (index_type)
            case (INDICES_TWOCENTER)
              if (m1tmp /= m2tmp) cycle
            case (INDICES_TWOCENTER_DIPX)
              if (abs(abs(m1tmp) - abs(m2tmp)) > 1 .or. &
                & (iand(indices_get_power_x(l1tmp, m1tmp) + indices_get_power_x(l2tmp, m2tmp), 1) == 0) .or. &
                & (iand(indices_get_power_y(l1tmp, m1tmp) + indices_get_power_y(l2tmp, m2tmp), 1) == 1)) cycle
            case (INDICES_TWOCENTER_DIPY)
              if (abs(abs(m1tmp) - abs(m2tmp)) > 1 .or. &
                & (iand(indices_get_power_x(l1tmp, m1tmp) + indices_get_power_x(l2tmp, m2tmp), 1) == 1) .or. &
                & (iand(indices_get_power_y(l1tmp, m1tmp) + indices_get_power_y(l2tmp, m2tmp), 1) == 0)) cycle
            case (INDICES_TWOCENTER_COULOMB)
              if (m1tmp /= 0 .or. m2tmp /= 0) cycle
            case default
              cycle
            end select
            index_max = index_max + 1
            buffer(1, index_max) = issh
            buffer(2, index_max) = jssh
            buffer(3, index_max) = l1tmp
            buffer(4, index_max) = l2tmp
            buffer(5, index_max) = m1tmp
            buffer(6, index_max) = m2tmp
          end do
        end do
      end do
    end do
    if (index_max == 0) return
    allocate (s1(index_max), s2(index_max), l1(index_max), l2(index_max), m1(index_max), m2(index_max))
    do ix = 1, index_max
      s1(ix) = buffer(1, ix)
      s2(ix) = buffer(2, ix)
      l1(ix) = buffer(3, ix)
      l2(ix) = buffer(4, ix)
      m1(ix) = buffer(5, ix)
      m2(ix) = buffer(6, ix)
    end do
    if (present(names)) then
      allocate (names(index_max))
      do ix = 1, index_max
        auxname1 = indices_get_name(l1(ix), m1(ix))
        auxname2 = indices_get_name(l2(ix), m2(ix))
        names(ix) = trim(auxname1)//trim(auxname2)
      end do
    end if
  end subroutine indices_twocenter_set

  pure integer function indices_get_power_x(l, m)
    integer, intent(in) :: l, m
    select case (l)
    case (0)
      indices_get_power_x = 0
    case (1)
      select case (m)
      case (1)
        indices_get_power_x = 1
      case default
        indices_get_power_x = 0
      end select
    case (2)
      select case (m)
      case (2)
        indices_get_power_x = 2
      case (-2, 1)
        indices_get_power_x = 1
      case default
        indices_get_power_x = 0
      end select
    case (3)
      select case (m)
      case (3)
        indices_get_power_x = 3
      case (-3, 2)
        indices_get_power_x = 2
      case (-2, 1)
        indices_get_power_x = 1
      case default
        indices_get_power_x = 0
      end select
    case default
      indices_get_power_x = 0
    end select
  end function indices_get_power_x

  pure integer function indices_get_power_y(l, m)
    integer, intent(in) :: l, m
    select case (l)
    case (0)
      indices_get_power_y = 0
    case (1)
      select case (m)
      case (-1)
        indices_get_power_y = 1
      case default
        indices_get_power_y = 0
      end select
    case (2)
      select case (m)
      case (2)
        indices_get_power_y = 2
      case (-2, -1)
        indices_get_power_y = 1
      case default
        indices_get_power_y = 0
      end select
    case (3)
      select case (m)
      case (-3)
        indices_get_power_y = 3
      case (2, 3)
        indices_get_power_y = 2
      case (-2, -1)
        indices_get_power_y = 1
      case default
        indices_get_power_y = 0
      end select
    case default
      indices_get_power_y = 0
    end select
  end function indices_get_power_y

  pure character(32) function indices_get_name(l, m, omit_suffix)
    integer, intent(in) :: l, m
    logical, intent(in), optional :: omit_suffix
    if (m < -l .or. m > l) then
      indices_get_name = ""
      return
    end if
    select case (l)
    case (0)
      indices_get_name = snames(m)
    case (1)
      indices_get_name = pnames(m)
    case (2)
      indices_get_name = dnames(m)
    case (3)
      indices_get_name = fnames(m)
    case default
      indices_get_name = " "
    end select
    if (present(omit_suffix)) then
      if (omit_suffix) indices_get_name = indices_get_name(1:1)
    end if
  end function indices_get_name

end module indices
