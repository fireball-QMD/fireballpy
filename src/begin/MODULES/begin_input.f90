module begin_input
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
integer :: ioption = 9  ! Will not be set from python and should not be used (hopefully)
integer :: nssh, nexcite, nznuc, nzval, nzval_ion, ioptim
real(kind=dp) :: nzval_pp
character(len=6) :: ppfile
character(len=8) :: ppionfile
character(len=1000) :: outpath
logical, dimension(:), allocatable :: sav
integer, dimension(:), allocatable :: lam
real(kind=dp), dimension(:), allocatable :: rcutoff, xocc, xocc0, xocc_ion, cmix, r0, v0, exc_weights
character(32), dimension(:), allocatable :: exc_names
end module begin_input
