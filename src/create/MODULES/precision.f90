module precision
  use iso_fortran_env, only: sp => real32, dp => real64
  implicit none(type, external)
  public
  integer, parameter :: wp = dp
end module precision
