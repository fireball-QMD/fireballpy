subroutine generate_basis (verbose_in)
  use constants, only: verbose
  implicit none
  logical, intent(in) :: verbose_in
  character(len=70) :: signature

  verbose = verbose_in
  signature = 'FireballPy'
  call create(signature)
end subroutine
