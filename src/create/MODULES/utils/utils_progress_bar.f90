submodule (utils) utils_progress_bar
  use iso_fortran_env, only: stderr => error_unit
  use precision, only: wp
  implicit none

  integer :: max_iter, bar_len, extra_length
  real(kind=wp) :: time
  character(32) :: fmt

contains

  module subroutine utils_progress_bar_prepare(max_iterations, length)
    integer, intent(in) :: max_iterations, length
    integer :: div, val, digits
    real(kind=wp) :: t
    character(4) :: cdigits
    digits = 0
    div = 1
    val = max_iterations
    do while (val > 0)
      val = val / div
      div = div*10
      digits = digits + 1
    end do
    extra_length = 35 + 2*digits
    bar_len = length - extra_length
    if (bar_len < 1) then
      write (stderr, '(a,i4,a,i4,a)')                                                &
      &  '[WARNING]: utils_progress_bar_prepare recieved too small length (got ',    &
      &  length, 'expected > ', extra_length, ')'
      return
    end if
    write (cdigits, '(i4)') digits
    fmt = '(a,i'//trim(cdigits)//',a,i'//trim(cdigits)//')'
    max_iter = max_iterations
    call cpu_time(time)
  end subroutine utils_progress_bar_prepare

  module subroutine utils_progress_bar_print(iteration, unit)
    integer, intent(in) :: iteration, unit
    integer :: i, percentage, len_done, mins, secs, minsr, secsr
    real(kind=wp) :: itps, newtime
    if (bar_len < 1) return
    percentage = (100*iteration)/max_iter
    len_done = nint(bar_len*percentage/100.0)
    call cpu_time(newtime)
    if (iteration > 1) then
      do i = 1, bar_len + extra_length
        write (unit, '(a)', advance='no') achar(8)
      end do
      itps = (iteration - 1)/(newtime - time)
      secs = nint(newtime - time)
      mins = secs/60
      secs = secs - 60*mins
      secsr = nint((max_iter - iteration)/itps)
      minsr = secsr/60
      secsr = secsr - 60*minsr
    else
      secs = 0
      secs = 0
      itps = 0.0_wp
      mins = 0
      secs = 0
      secsr = 0
      minsr = 0
      secsr = 0
    end if
    write (unit, '(i3,a)', advance='no') percentage, '% |'
    do i = 1, len_done
      write (unit, '(a)', advance='no') '#'
    end do
    do i = len_done + 1, bar_len
      write (unit, '(a)', advance='no') ' '
    end do
    write (unit, fmt, advance='no') '| ', iteration, '/', max_iter
    write (unit, '(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,f6.2,a)', advance='no')       &
    &  ' [', mins, ':', secs, '<', minsr, ':', secsr, ', ', itps, 'it/s]'
    flush(unit)
  end subroutine utils_progress_bar_print

  module subroutine utils_progress_bar_clear(unit)
    integer, intent(in) :: unit
    integer :: i
    do i = 1, bar_len + extra_length
      write (unit, '(a)', advance='no') achar(8)
      write (unit, '(a)', advance='no') ' '
      write (unit, '(a)', advance='no') achar(8)
    end do
    flush(unit)
  end subroutine utils_progress_bar_clear

end submodule utils_progress_bar
