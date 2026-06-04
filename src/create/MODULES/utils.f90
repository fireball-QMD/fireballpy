module utils
  implicit none
  private
  public :: utils_progress_bar_prepare, utils_progress_bar_print, utils_progress_bar_clear

  interface
    module subroutine utils_progress_bar_prepare(max_iterations, length)
      implicit none
      integer, intent(in) :: max_iterations, length
    end subroutine utils_progress_bar_prepare
  end interface

  interface
    module subroutine utils_progress_bar_print(iteration, unit)
      implicit none
      integer, intent(in) :: iteration, unit
    end subroutine utils_progress_bar_print
  end interface

  interface
    module subroutine utils_progress_bar_clear(unit)
      implicit none
      integer, intent(in) :: unit
    end subroutine utils_progress_bar_clear
  end interface
end module utils
