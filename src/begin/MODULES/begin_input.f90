module begin_input
  use precision
  implicit none
  integer :: ioption, nexcite, nssh, nznuc, nzval, nzval_ion, ioptim
  real(kind=long) :: nzval_pp
  character(len=10) :: atomname
  character(len=6) :: ppfile
  character(len=8) :: ppionfile
  character(len=1000) :: outpath
  logical, dimension(:), allocatable :: sav
  integer, dimension(:), allocatable :: lam
  real(kind=long), dimension(:), allocatable :: a0, rcutoff, xocc, xocc0, xocc_ion, cmix, r0, v0
  character(len=11), dimension(:), allocatable :: filename_wf
  character(len=11), dimension(:), allocatable :: filename_na
  character(len=12), dimension(:), allocatable :: filename_ewf, filename_ena
end module begin_input
