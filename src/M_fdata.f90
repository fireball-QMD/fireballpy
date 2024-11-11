module M_fdata
  use iso_c_binding

  !====================================
  ! THE CODE SHOULD RUN EVENTUALLY WITHOUT THESE VARIABLES
  integer(c_long) :: itheory = 1 !DOGS 
  integer(c_long) :: itheory_xc = 2 !McWEDA
  logical :: debug = .True.
  character(len=100) :: infofname = 'info.py.dat'
  !===================================
  
  ! info.dat variables
  integer(c_long) :: nsh_max, nshPP_max, isorpmax, isorpmax_xc, nspecies
  character (len=1000) fdataLocation
  integer(c_long), dimension (:), allocatable :: nzx, nssh, nsshPP
  real(c_double), dimension (:), allocatable :: etotatom, smass, rc_PP
  character (len=2), dimension (:), allocatable :: symbolA
  integer(c_long), dimension (:,:), allocatable :: lssh, lsshPP
  real(c_double), dimension (:,:), allocatable :: rcutoff, cl_PP, Qneutral
  character (len=25), dimension (:,:), allocatable :: wavefxn, napot

  ! Maximum number of two-center matrix elements: (calculated in make_munu.f90)
  ! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 39
  integer(c_long) :: ME2c_max, ME2cPP_max, ME2cDipY_max, ME2cDipX_max

  ! Maximum number of three-center matrix elements: (calculated in make_munu.f90)
  ! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
  integer(c_long) :: ME3c_max

  ! Maximum number of two and three-center matrix elements in spherical density
  ! approximation (OLSXC) (calculated in make_munuS.f90)
  ! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
  integer(c_long) :: MES_max

  integer(c_long), dimension (:), allocatable :: num_orb
  integer(c_long), dimension (:,:), allocatable :: index_max2c, index_max2cDipX, index_max2cDipY, index_max3c
  integer(c_long), dimension (:,:,:), allocatable :: mu, nu, mvalue, muDipX, nuDipX, muDipY, nuDipY

  ! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
  integer(c_long), dimension (:), allocatable :: num_orbPP
  integer(c_long), dimension (:,:), allocatable :: index_maxPP
  integer(c_long), dimension (:,:,:), allocatable :: muPP, nuPP

  ! These variables are specifically for spherical density approximation 
  ! used in OLSXC method
  integer(c_long), dimension (:,:), allocatable :: index_maxS
  integer(c_long), dimension (:,:,:), allocatable :: muS, nuS, mvalueS

  ! new Intra-atomic Dipole: One-center case (for the time being)
  integer(c_long), dimension(:,:), allocatable :: muR, nuR, alphaR, betaR
  real(c_double), dimension(:,:), allocatable :: IR

  !TODO idipole, icluster, siempre lee interaccion 10 y 11

  ! One center integrals
  character (len=9), dimension (3), parameter :: onecfname = (/'xc1c_dqi ','nuxc1crho','exc1crho '/)
  real(c_double), dimension (:,:), allocatable :: exc1c_0
  real(c_double), dimension (:,:,:), allocatable :: xcnu1c
  real(c_double), dimension (:,:,:), allocatable :: xcnu1cs, exc1c0, nuxc1c, d2exc1c
  real(c_double), dimension (:,:,:,:), allocatable :: exc1c, dexc1c, dnuxc1c
  real(c_double), dimension (:,:,:,:,:), allocatable :: d2nuxc1c

  ! Two center integrals
  integer(c_long), parameter :: nfofx = 207 ! AQUI
  character (len=11), dimension (23), parameter :: twocfname = (/'overlap    ','vna_ontopl ','vna_ontopr ', &
    & 'vna_atom   ','vnl        ','xc_ontop   ','xc_atom    ','xc_corr    ','dipole_z   ','dipole_y   ','dipole_x   ', &
    & 'coulomb    ','kinetic    ','nuxc       ','den_ontopl ','den_ontopr ','den_atom   ','dnuxc_ol   ','dnuxc_or   ', &
    & 'denS_ontopl','denS_ontopr','denS_atom  ','overlapS   '/)
  integer(c_long) :: interactions2c_max
  integer(c_long), dimension (1:23,0:8) :: ind2c
  integer(c_long), dimension (:,:,:), allocatable :: numz2c
  real(c_double), dimension (:,:,:), allocatable :: z2cmax
  real(c_double), dimension (:,:,:,:,:,:), allocatable :: splineint_2c

  ! Three-center integrals
  integer(c_long), parameter :: numXmax = 31 ! AQUI
  integer(c_long), parameter :: numYmax = 31 ! AQUI
  integer(c_long), parameter :: ntheta = 5
  character (len=4), dimension (4), parameter :: threecfname = (/'bcna','xc3c','den3','deS3'/)
  !integer(c_long), parameter :: ideriv_max = 6
  integer(c_long), dimension (:,:,:), allocatable :: icon3c
  ! Neutral (charged) atom interactions; there are five bcna arrays - one for each theta
  integer(c_long), dimension (:,:), allocatable :: numx3c_bcna, numy3c_bcna
  real(c_double), dimension (:,:), allocatable :: hx_bcna, hy_bcna, x3cmax_bcna, y3cmax_bcna
  real(c_double), dimension (:,:,:,:,:), allocatable :: bcna_01, bcna_02, bcna_03, bcna_04, bcna_05
  ! XC interactions; 7 implies different derivative types; there are five xc3c arrays - one for each theta
  !integer(c_long), dimension (:,:), allocatable :: numx3c_xc3c, numy3c_xc3c
  !real(c_double), dimension (:,:), allocatable :: hx_xc3c, hy_xc3c, x3cmax_xc3c, y3cmax_xc3c
  !real(c_double), dimension (:,:,:,:,:), allocatable :: xc3c_01, xc3c_02, xc3c_03, xc3c_04, xc3c_05
  ! XC interactions; 7 implies different derivative types; there are five den3 arrays - one for each theta, used only for SNXC method
  integer(c_long), dimension (:,:), allocatable :: numx3c_den3, numy3c_den3
  real(c_double), dimension (:,:), allocatable :: hx_den3, hy_den3, x3cmax_den3, y3cmax_den3
  real(c_double), dimension (:,:,:,:,:), allocatable :: den3_01, den3_02, den3_03, den3_04, den3_05
  real(c_double), dimension (:,:,:,:,:), allocatable :: den3S_01, den3S_02, den3S_03, den3S_04, den3S_05

contains
    function append_string (filename, extension)
        implicit none
        character (len = 200) append_string
        character (len = 200) filename
        character (len = 200) extension
        integer(c_long) lenx
        integer(c_long) leny
        integer(c_long) len_trim
        lenx = len_trim(filename)
        if (lenx .eq. 0) then
            append_string = ' '
        else
            leny = len_trim(extension)
            if (leny .eq. 0) then
                append_string = filename(1:lenx)
            else
                append_string = filename(1:lenx)//extension(1:leny)
            end if
        end if
        return
    end function append_string

end module M_fdata
