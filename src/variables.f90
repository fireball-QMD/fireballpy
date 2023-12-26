module variables
  real :: num1 = 0
  real :: num2 = 0
  real :: resultado = 0

  integer, parameter :: norbitals = 4 
  complex*16 :: yyyy(norbitals,norbitals), eigen(norbitals)
  complex*16, allocatable, dimension (:) :: work
  real*8, allocatable, dimension (:) :: rwork
  integer lwork, lrwork
  integer :: i, j
  character(len=20) :: filename_in = 'fireball.in'

end module variables

