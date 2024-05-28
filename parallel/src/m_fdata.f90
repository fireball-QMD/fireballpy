module M_fdata

  !====================================
  ! THE CODE SHOULD RUN EVENTUALLY WITHOUT THESE VARIABLES
  integer :: itheory = 1 !DOGS 
  integer :: itheory_xc = 2 !McWEDA
  logical :: debug = .True.
  integer :: fdata_is_load = 0
  !===================================

  !===========================================
  ! VARIABLES WHICH WILL BE PASSED FROM PYTHON
  !f2py start inblock
  integer :: nsh_max, nspecies, ME2c_max, ME2cPP_max, ME2cDipY_max, ME2cDipX_max, ME3c_max, MES_max
  character (len=1000) :: fdataLocation
  integer, dimension(nspecies) :: nzx
  integer, dimension (nsh_max) :: nssh, nsshPP
  real*8, dimension (nspecies) :: etotatom, smass, rc_PP
  character (len=2), dimension (nspecies) :: symbolA
  integer, dimension (nsh_max,nspecies) :: lssh, lsshPP
  real*8, dimension (nspecies,nsh_max) :: rcutoff
  real*8, dimension (nsh_max,nspecies) :: Qneutral
  character (len=25), dimension (nsh_max,nspecies) :: wavefxn
  character (len=25), dimension (0:nsh_max,nspecies) :: napot
  !f2py end inblock
  !===========================================

  ! Maximum number of two-center matrix elements: (calculated in make_munu.f90)
  ! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 39
  !integer :: ME2c_max, ME2cPP_max, ME2cDipY_max, ME2cDipX_max

  ! Maximum number of three-center matrix elements: (calculated in make_munu.f90)
  ! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
  !integer :: ME3c_max

  ! Maximum number of two and three-center matrix elements in spherical density
  ! approximation (OLSXC) (calculated in make_munuS.f90)
  ! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
  !integer :: MES_max

  integer, dimension (nspecies) :: num_orb
  integer, dimension (nspecies,nspecies) :: index_max2c, index_max2cDipX, index_max2cDipY, index_max3c
  integer, dimension (ME3c_max,nspecies,nspecies) :: mu, nu, mvalue
  integer, dimension (ME2cDipX_max,nspecies,nspecies) :: muDipX, nuDipX
  integer, dimension (ME2cDipY_max,nspecies,nspecies) :: muDipY, nuDipY

  ! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
  integer, dimension (nspecies) :: num_orbPP
  integer, dimension (nspecies,nspecies) :: index_maxPP
  integer, dimension (ME2cPP_max,nspecies,nspecies) :: muPP, nuPP

  ! These variables are specifically for spherical density approximation 
  ! used in OLSXC method
  integer, dimension (nspecies,nspecies) :: index_maxS
  integer, dimension (MES_max,nspecies,nspecies) :: muS, nuS, mvalueS

  ! new Intra-atomic Dipole: One-center case (for the time being)
  !integer, dimension(:,:), allocatable :: muR, nuR, alphaR, betaR
  !real*8, dimension(:,:), allocatable :: IR

  !TODO idipole, icluster, siempre lee interaccion 10 y 11

  ! One center integrals
  character (len=9), dimension (3) :: onecfname
  real*8, dimension (nspecies,nsh_max,nsh_max) :: exc1c0, nuxc1c
  real*8, dimension (nspecies,nsh_max,nsh_max,nsh_max) :: dexc1c, dnuxc1c
  !real*8, dimension (nspecies,nsh_max,nsh_max) :: d2exc1c
  !real*8, dimension (nspecies,nsh_max,nsh_max,nsh_max,nsh_max) :: d2nuxc1c
  !real*8, dimension (:,:), allocatable :: exc1c_0
  !real*8, dimension (:,:,:), allocatable :: xcnu1c, xcnu1cs
  !real*8, dimension (:,:,:,:), allocatable :: exc1c

  ! Two center integrals
  integer, parameter :: nfofx = 207 ! AQUI
  character (len=11), dimension (23) :: twocfname
  integer :: errno2c
  integer :: interactions2c_max
  integer, dimension (1:23,0:8) :: ind2c
  integer, dimension (27+11*nsh_max,nspecies,nspecies) :: numz2c
  real*8, dimension (nspecies,nsh_max) :: cl_PP
  real*8, dimension (27+11*nsh_max,nspecies,nspecies) :: z2cmax
  real*8, dimension (4,ME2c_max,nfofx,27+11*nsh_max,nspecies,nspecies) :: splineint_2c

  ! Three-center integrals
  integer, parameter :: numXmax = 31 ! AQUI
  integer, parameter :: numYmax = 31 ! AQUI
  integer, parameter :: ntheta = 5
  character (len=4), dimension (4) :: threecfname
  !integer, parameter :: ideriv_max = 6
  integer :: errno3c
  integer, dimension (nspecies,nspecies,nspecies) :: icon3c
  ! Neutral (charged) atom interactions; there are five bcna arrays - one for each theta
  integer, dimension (0:nsh_max,nspecies*nspecies*nspecies) :: numx3c_bcna, numy3c_bcna
  real*8, dimension (0:nsh_max,nspecies*nspecies*nspecies) :: hx_bcna, hy_bcna, x3cmax_bcna, y3cmax_bcna
  real*8, dimension (numXmax,numYmax,ME3c_max,0:nsh_max,nspecies*nspecies*nspecies) :: bcna_01, bcna_02, bcna_03, bcna_04, bcna_05
  ! XC interactions; 7 implies different derivative types; there are five xc3c arrays - one for each theta
  !integer, dimension (0:ideriv_max,nspecies*nspecies*nspecies), allocatable :: numx3c_xc3c, numy3c_xc3c
  !real*8, dimension (0:ideriv_max,nspecies*nspecies*nspecies), allocatable :: hx_xc3c, hy_xc3c, x3cmax_xc3c, y3cmax_xc3c
  !real*8, dimension (numXmax,numYmax,ME3c_max,0:ideriv_max,nspecies*nspecies*nspecies) :: xc3c_01, xc3c_02, xc3c_03, xc3c_04, xc3c_05
  ! XC interactions; 7 implies different derivative types; there are five den3 arrays - one for each theta, used only for SNXC method
  integer, dimension (0:nsh_max,nspecies*nspecies*nspecies) :: numx3c_den3, numy3c_den3
  real*8, dimension (0:nsh_max,nspecies*nspecies*nspecies) :: hx_den3, hy_den3, x3cmax_den3, y3cmax_den3
  real*8, dimension (numXmax,numYmax,ME3c_max,0:nsh_max,nspecies*nspecies*nspecies) :: den3_01, den3_02, den3_03, den3_04, den3_05
  real*8, dimension (numXmax,numYmax,ME3c_max,0:nsh_max,nspecies*nspecies*nspecies) :: den3S_01, den3S_02, den3S_03, den3S_04, den3S_05
end module M_fdata
