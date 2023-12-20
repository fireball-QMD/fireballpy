module M_fdata

  logical :: debug = .True.
  integer, parameter :: nsh_max = 6
  integer :: nspecies
  character (len = 200) fdataLocation
  integer, dimension (:), allocatable :: nzx
  character (len = 2), dimension (:), allocatable :: symbolA
  real, dimension (:), allocatable :: etotatom
  real, dimension (:), allocatable :: smass
  real, dimension (:), allocatable :: rc_PP
  real, dimension (:,:), allocatable :: rcutoff
  real, dimension (:,:), allocatable :: cl_PP
  integer, dimension (:), allocatable :: nssh
  integer, dimension (:), allocatable :: nsshPP
  integer, dimension (:,:), allocatable :: lssh
  integer, dimension (:,:), allocatable :: lsshPP
  real, dimension (:,:), allocatable :: Qneutral
  character (len=25), dimension (:,:), allocatable :: wavefxn
  character (len=25), dimension (:,:), allocatable :: napot
   

end module M_fdata
