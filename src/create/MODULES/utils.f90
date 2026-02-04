module utils
  implicit none
  public

contains

  subroutine utils_progress_bar(iteration, max_iterations, length, unit)
    integer, intent(in) :: iteration, max_iterations, length, unit
    integer :: i, percentage, len_done
    percentage = (100*iteration)/max_iterations
    len_done = nint(length*percentage/100.0)
    if (iteration > 1) then
      do i = 1, length + 17
        write (unit, '(a)', advance='no') achar(8)
      end do
    end if
    write (unit, '(i3,a)', advance='no') percentage, '% |'
    do i = 1, len_done
      write (unit, '(a)', advance='no') '#'
    end do
    do i = len_done + 1, length
      write (unit, '(a)', advance='no') ' '
    end do
    write (unit, '(a,i4,a,i4)', advance='no') '| ', iteration, '/', max_iterations
    flush(unit)
  end subroutine utils_progress_bar

  subroutine utils_clean_progress_bar(length, unit)
    integer, intent(in) :: length, unit
    integer :: i
    do i = 1, length + 17
      write (unit, '(a)', advance='no') achar(8)
    end do
    do i = 1, length + 17
      write (unit, '(a)', advance='no') ' '
    end do
    do i = 1, length + 17
      write (unit, '(a)', advance='no') achar(8)
    end do
    flush(unit)
  end subroutine utils_clean_progress_bar

end module utils
