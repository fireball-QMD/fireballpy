module utils
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use :: constants, only:err_len
  implicit none
  private
  public :: utils_open

contains

  integer function utils_open(filename, mode)
    character(*), intent(in) :: filename
    character(*), intent(in), optional :: mode
    integer :: stat
    character(3) :: mode_
    character(err_len) :: errmsg
    character(:), allocatable :: action_, position_, status_, access_, form_

    if (.not. present(mode)) then
      mode_ = "r t"
    else
      select case (len(mode))
      case (1)
        mode_ = mode//" t"
      case (2)
        mode_ = mode//"t"
      case default
        mode_ = trim(adjustl(mode(1:3)))
      end select
    end if

    select case (mode_(1:2))
    case ("r")
      action_ = "read"
      position_ = "asis"
      status_ = "old"
    case ("w")
      action_ = "write"
      position_ = "asis"
      status_ = "replace"
    case ("a")
      action_ = "write"
      position_ = "append"
      status_ = "old"
    case ("x")
      action_ = "write"
      position_ = "asis"
      status_ = "new"
    case ("r+")
      action_ = "readwrite"
      position_ = "asis"
      status_ = "old"
    case ("w+")
      action_ = "readwrite"
      position_ = "asis"
      status_ = "replace"
    case ("a+")
      action_ = "readwrite"
      position_ = "append"
      status_ = "old"
    case ("x+")
      action_ = "readwrite"
      position_ = "asis"
      status_ = "new"
    case default
      errmsg = "[ERROR] Unsupported mode: "//mode_(1:2)
      write (stderr, "(a)") trim(errmsg)
      error stop 1
    end select

    select case (mode_(3:3))
    case ("t")
      form_ = "formatted"
      access_ = "sequential"
    case ("b")
      form_ = "unformatted"
      access_ = "stream"
    case default
      errmsg = "[ERROR] Unsupported mode: "//mode_(3:3)
      write (stderr, "(a)") trim(errmsg)
      error stop 1
    end select

    open (newunit=utils_open, file=trim(filename), action=action_, position=position_, status=status_, &
      &   access=access_, form=form_, iostat=stat, iomsg=errmsg)
    if (stat /= 0) then
      write (stderr, "(a)") "[ERROR] Could not open "//trim(filename)//": "//trim(errmsg)
      error stop 1
    end if
  end function utils_open

end module utils
