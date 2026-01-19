submodule (math) math_lstsq_impl
  use iso_fortran_env, only: sp => real32, dp => real64
  implicit none

contains

  module subroutine slstsq(a, b, is_a_trans, info)
    real(kind=sp), intent(in out) :: a(:,:), b(:,:)
    logical, intent(in), optional :: is_a_trans
    integer, intent(out), optional :: info
    integer :: nrows, ncols, nrhs, nobs, lwork, stat
    character(1) :: ctrans
    real(kind=sp), allocatable :: work(:)

    nrows = size(a, 1)
    ncols = size(a, 2)
    nobs = size(b, 1)
    nrhs =  size(b, 2)
    if (.not. present(is_a_trans)) then
      ctrans = 'N'
    else if (is_a_trans) then
      ctrans = 'T'
    else
      ctrans = 'N'
    end if
    if (((ctrans == 'N') .and. (nrows /= nobs)) .or. ((ctrans == 'T') .and. (ncols /= nobs))) then
      if (present(info)) info = -1
      return
    end if

    lwork = min(nrows, ncols) + max(min(nrows, ncols), nrhs)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate(work(lwork), stat=stat)
    if(present(info)) info = stat
    if(stat /= 0) return
    call sgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, stat)
    deallocate(work)
    if(present(info)) info = stat
    if(stat /= 0) return
  end subroutine slstsq

  module subroutine dlstsq(a, b, is_a_trans, info)
    real(kind=dp), intent(in out) :: a(:,:), b(:,:)
    logical, intent(in), optional :: is_a_trans
    integer, intent(out), optional :: info
    integer :: nrows, ncols, nrhs, nobs, lwork, stat
    character(1) :: ctrans
    real(kind=dp), allocatable :: work(:)

    nrows = size(a, 1)
    ncols = size(a, 2)
    nobs = size(b, 1)
    nrhs =  size(b, 2)
    if (.not. present(is_a_trans)) then
      ctrans = 'N'
    else if (is_a_trans) then
      ctrans = 'T'
    else
      ctrans = 'N'
    end if
    if (((ctrans == 'N') .and. (nrows /= nobs)) .or. ((ctrans == 'T') .and. (ncols /= nobs))) then
      if (present(info)) info = -1
      return
    end if

    lwork = min(nrows, ncols) + max(min(nrows, ncols), nrhs)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate(work(lwork), stat=stat)
    if(present(info)) info = stat
    if(stat /= 0) return
    call dgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, stat)
    deallocate(work)
    if(present(info)) info = stat
    if(stat /= 0) return
  end subroutine dlstsq
end submodule math_lstsq_impl
