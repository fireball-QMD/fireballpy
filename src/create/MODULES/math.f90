module math
  use iso_fortran_env, only: sp => real32, dp => real64
  implicit none
  private

  public :: math_lstsq
  public :: math_Ylm

  interface math_lstsq
    module integer function slstsq(a, b, is_a_trans)
      implicit none
      real(kind=sp), intent(in out) :: a(:,:), b(:,:)
      logical, intent(in), optional :: is_a_trans
    end function slstsq
    module integer function dlstsq(a, b, is_a_trans)
      implicit none
      real(kind=dp), intent(in out) :: a(:,:), b(:,:)
      logical, intent(in), optional :: is_a_trans
    end function dlstsq
  end interface math_lstsq

  interface math_Ylm
    module real(kind=sp) function sYlm(l, m, prec)
      implicit none
      integer, intent(in) :: l, m
      real(kind=sp), intent(in) :: prec
    end function sYlm
    module real(kind=dp) function dYlm(l, m, prec)
      implicit none
      integer, intent(in) :: l, m
      real(kind=dp), intent(in) :: prec
    end function dYlm
  end interface math_Ylm

end module math
