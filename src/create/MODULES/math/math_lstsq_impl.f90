submodule (math) math_lstsq_impl
  implicit none

contains

  module integer function slstsq(a, b, is_a_trans)
    real(kind=sp), intent(in out) :: a(:,:), b(:,:)
    logical, intent(in), optional :: is_a_trans
    integer :: nrows, ncols, nrhs, nobs, lwork
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
      slstsq = -1
      return
    end if

    lwork = min(nrows, ncols) + max(min(nrows, ncols), nrhs)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate(work(lwork), stat=slstsq)
    if(slstsq /= 0) return
    call sgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, slstsq)
    deallocate(work)
    if(slstsq /= 0) return
    slstsq = 0
    return
  end function slstsq

  module integer function dlstsq(a, b, is_a_trans)
    real(kind=dp), intent(in out) :: a(:,:), b(:,:)
    logical, intent(in), optional :: is_a_trans
    integer :: nrows, ncols, nrhs, nobs, lwork
    character(1) :: ctrans
    real(kind=dp), allocatable :: work(:)

    nrows = size(a, 1)
    ncols = size(a, 2)
    nobs = size(b, 1)
    nrhs = size(b, 2)
    if (.not. present(is_a_trans)) then
      ctrans = 'N'
    else if (is_a_trans) then
      ctrans = 'T'
    else
      ctrans = 'N'
    end if
    if (((ctrans == 'N') .and. (nrows /= nobs)) .or. ((ctrans == 'T') .and. (ncols /= nobs))) then
      dlstsq = -1
      return
    end if

    lwork = min(nrows, ncols) + max(min(nrows, ncols), nrhs)
    lwork = ibset(0, bit_size(lwork) - leadz(lwork))
    allocate(work(lwork), stat=dlstsq)
    if(dlstsq /= 0) return
    call dgels(ctrans, nrows, ncols, nrhs, a, nrows, b, nobs, work, lwork, dlstsq)
    deallocate(work)
    if(dlstsq /= 0) return
    dlstsq = 0
    return
  end function dlstsq

end submodule math_lstsq_impl
