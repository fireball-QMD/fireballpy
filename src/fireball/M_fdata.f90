module M_fdata
  use, intrinsic :: iso_fortran_env, only: double => real64
  use math, only: math_interp_t
  implicit none
  public

  ! Where the FData tables live. The rest of what used to come out of
  ! 'info.dat' is now set from Python through 'set_fdata_info'.
  character(len=400) :: fdataLocation

  integer :: nsh_max, nshPP_max, isorpmax, isorpmax_xc, nspecies

  integer, dimension (:), allocatable :: nzx, nssh, nsshPP
  real(double), dimension (:), allocatable :: rc_PP
  integer, dimension (:,:), allocatable :: lssh, lsshPP
  real(double), dimension (:,:), allocatable :: rcutoff, cl_PP, Qneutral, Qref

  ! Maximum number of two-center matrix elements
  ! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 39
  integer :: ME2c_max, ME2cPP_max, ME2cDipY_max, ME2cDipX_max

  ! Maximum number of three-center matrix elements
  ! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
  integer :: ME3c_max

  ! Maximum number of two and three-center matrix elements in spherical density approximation (OLSXC)
  ! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
  integer :: MES_max

  integer, dimension(:), allocatable :: num_orb
  integer :: numorb_max
  integer, dimension(:, :), allocatable :: index_max2c, index_max2cDipX, index_max2cDipY, index_max3c
  integer, dimension(:, :, :), allocatable :: mu, nu, mvalue, muDipX, nuDipX, muDipY, nuDipY

  ! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
  integer, dimension(:), allocatable :: num_orbPP
  integer, dimension(:, :), allocatable :: index_maxPP
  integer, dimension(:, :, :), allocatable :: muPP, nuPP

  ! These variables are specifically for spherical density approximation
  ! used in OLSXC method
  integer, dimension(:, :), allocatable :: index_maxS
  integer, dimension(:, :, :), allocatable :: muS, nuS, mvalueS

  character(len=32), allocatable :: exc_names(:)
  real(double), allocatable :: exc_weights(:)

  !TODO idipole, icluster, siempre lee interaccion 10 y 11

  ! One center integrals
  character (len=2), dimension (1), parameter :: onecfname = ["xc"]
  real(double), dimension (:,:), allocatable :: exc_1c_0
  real(double), dimension (:,:,:), allocatable :: vxc_1c_0, fxc_1c
  real(double), dimension (:,:,:,:), allocatable :: gxc_1c


  character(len=11), dimension(20), parameter :: twocfname = [ &
                                                 "overlap    ", "vna_ontopl ", "vna_ontopr ", "vna_atom   ", "vnl        ", &
                                                 "vxc_refchs ", "vxc_ontopl ", "vxc_ontopr ", "dipole_z   ", "dipole_y   ", &
                                                 "dipole_x   ", "coulomb    ", "kinetic    ", "den_ontopl ", "den_ontopr ", &
                                                 "den_atom   ", "denS_ontopl", "denS_ontopr", "denS_atom  ", "overlapS   "]
  integer, parameter :: TWOCENTER_OVERLAP = 1
  integer, parameter :: TWOCENTER_VNA_L = 2
  integer, parameter :: TWOCENTER_VNA_R = 3
  integer, parameter :: TWOCENTER_VNA_A = 4
  integer, parameter :: TWOCENTER_VNL = 5
  integer, parameter :: TWOCENTER_VXC_0 = 6
  integer, parameter :: TWOCENTER_VXC_L = 7
  integer, parameter :: TWOCENTER_VXC_R = 8
  integer, parameter :: TWOCENTER_DIP_Z = 9
  integer, parameter :: TWOCENTER_DIP_Y = 10
  integer, parameter :: TWOCENTER_DIP_X = 11
  integer, parameter :: TWOCENTER_COULOMB = 12
  integer, parameter :: TWOCENTER_KINETIC = 13
  integer, parameter :: TWOCENTER_DEN_L = 14
  integer, parameter :: TWOCENTER_DEN_R = 15
  integer, parameter :: TWOCENTER_DEN_A = 16
  integer, parameter :: TWOCENTER_DENS_L = 17
  integer, parameter :: TWOCENTER_DENS_R = 18
  integer, parameter :: TWOCENTER_DENS_A = 19
  integer, parameter :: TWOCENTER_OVERLAPS = 20
  integer, parameter :: initype(20) = [ &
                        0, 0, 0, 0, 0, 0, & ! 1–6
                        1, 1, & ! 7–8
                        0, 0, 0, 0, 0, & ! 9–13
                        1, 1, 1, 1, 1, 1, & ! 14–19
                        0]                  ! 20
  integer, parameter :: maxtype_2c(20) = [ &
                        0, 1, 1, 1, 0, 0, & ! 1–6
                        1, 1, & ! 7–8
                        0, 0, 0, 0, 0, & ! 9–13
                        1, 1, 1, 1, 1, 1, & ! 14–19
                        0]                  ! 20

  integer, dimension(1:20, 0:8) :: ind2c
  type(math_interp_t), allocatable :: splineint_2c(:,:,:,:)

  ! Three-center integrals
  integer, parameter :: numXmax = 31 ! AQUI
  integer, parameter :: numYmax = 31 ! AQUI
  integer, parameter :: ntheta = 5
  character(len=4), dimension(4), parameter :: threecfname = ["bcna", "xc3c", "den3", "deS3"]
  !integer, parameter :: ideriv_max = 6
  integer, dimension(:, :, :), allocatable :: icon3c
  ! Neutral (charged) atom interactions; there are five bcna arrays - one for each theta
  integer, dimension(:, :), allocatable :: numx3c_bcna, numy3c_bcna
  real(double), dimension(:, :), allocatable :: hx_bcna, hy_bcna, x3cmax_bcna, y3cmax_bcna
  real(double), dimension(:, :, :, :, :), allocatable :: bcna_01, bcna_02, bcna_03, bcna_04, bcna_05
  integer, dimension(:, :), allocatable :: numx3c_den3, numy3c_den3
  real(double), dimension(:, :), allocatable :: hx_den3, hy_den3, x3cmax_den3, y3cmax_den3
  real(double), dimension(:, :, :, :, :), allocatable :: den3_01, den3_02, den3_03, den3_04, den3_05
  real(double), dimension(:, :, :, :, :), allocatable :: den3S_01, den3S_02, den3S_03, den3S_04, den3S_05

end module M_fdata
