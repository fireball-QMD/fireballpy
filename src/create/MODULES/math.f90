module math
  use iso_fortran_env, only: sp => real32, dp => real64
  implicit none(type, external)
  private

  public :: math_lstsq
  interface math_lstsq
    module subroutine slstsq(a, b, is_a_trans, info)
      implicit none(type, external)
      real(kind=sp), intent(in out) :: a(:,:), b(:,:)
      logical, intent(in), optional :: is_a_trans
      integer, intent(out), optional :: info
    end subroutine slstsq
    module subroutine dlstsq(a, b, is_a_trans, info)
      implicit none(type, external)
      real(kind=dp), intent(in out) :: a(:,:), b(:,:)
      logical, intent(in), optional :: is_a_trans
      integer, intent(out), optional :: info
    end subroutine dlstsq
  end interface math_lstsq
end module math
