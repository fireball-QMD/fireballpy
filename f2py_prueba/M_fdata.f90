module M_fdata

    use M_constants
    implicit none

    logical :: debug = .True. ! Tendremos q borrarla

    !integer :: itheory = 1 !DOGS 
    !integer :: itheory_xc = 2 !McWEDA
    !integer :: ispin = 0
    integer :: ideriv_max = 6

    ! load_fdata
    integer :: nsh_max
    integer :: nshPP_max
    integer :: isorpmax
    integer :: isorpmax_xc
    integer :: nspecies
    character (len = 200) fdataLocation

    integer, dimension (:), allocatable :: nzx
    integer, dimension (:), allocatable :: nssh
    integer, dimension (:), allocatable :: nsshPP
    real, dimension (:), allocatable :: etotatom
    real, dimension (:), allocatable :: smass
    real, dimension (:), allocatable :: rc_PP
    character (len = 2), dimension (:), allocatable :: symbolA

    integer, dimension (:,:), allocatable :: lssh
    integer, dimension (:,:), allocatable :: lsshPP
    real, dimension (:,:), allocatable :: rcutoff
    real, dimension (:,:), allocatable :: cl_PP
    real, dimension (:,:), allocatable :: Qneutral
    character (len=25), dimension (:,:), allocatable :: wavefxn
    character (len=25), dimension (:,:), allocatable :: napot

    ! Maximum number of two-center matrix elements: (calculated in make_munu.f90")
    ! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 39

    integer :: ME2c_max
    integer :: ME2cPP_max
    integer :: ME2cDipY_max
    integer :: ME2cDipX_max

    ! Maximum number of three-center matrix elements: (calculated in make_munu.f90")
    ! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
    integer :: ME3c_max

    ! Maximum number of two and three-center matrix elements in spherical density
    ! approximation (OLSXC) (calculated in make_munuS.f90")
    ! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
    integer :: MES_max

    integer, dimension (:, :), allocatable :: index_max2c
    integer, dimension (:, :), allocatable :: index_max2cDipY
    integer, dimension (:, :), allocatable :: index_max2cDipX
    integer, dimension (:, :), allocatable :: index_max3c

    integer, dimension (:, :, :), allocatable :: mu
    integer, dimension (:, :, :), allocatable :: nu
    integer, dimension (:, :, :), allocatable :: mvalue
    integer, dimension (:), allocatable :: num_orb

    integer, dimension (:, :, :), allocatable :: muDipY
    integer, dimension (:, :, :), allocatable :: nuDipY
    integer, dimension (:, :, :), allocatable :: muDipX
    integer, dimension (:, :, :), allocatable :: nuDipX


    ! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
    integer, dimension (:, :), allocatable :: index_maxPP
    integer, dimension (:, :, :), allocatable :: muPP
    integer, dimension (:, :, :), allocatable :: nuPP
    integer, dimension (:), allocatable :: num_orbPP


    ! These variables are specifically for spherical density approximation 
    ! used in OLSXC method
    integer, dimension (:, :), allocatable :: index_maxS
    integer, dimension (:, :, :), allocatable :: muS
    integer, dimension (:, :, :), allocatable :: nuS
    integer, dimension (:, :, :), allocatable :: mvalueS


    ! new Intra-atomic Dipole: One-center case (for the time being)
    integer, dimension(:,:), allocatable :: muR
    integer, dimension(:,:), allocatable :: nuR
    integer, dimension(:,:), allocatable :: alphaR
    integer, dimension(:,:), allocatable :: betaR
    real, dimension(:,:), allocatable :: IR
    integer, dimension(:), allocatable :: Nlines_vdip1c


    !TODO idipole, icluster, siempre lee interaccion 10 y 11

    ! Three-center integrals
    integer, dimension (:, :, :), allocatable :: icon3c

    ! Neutral (charged) atom interactions; there are five bcna arrays - one for each theta
    integer, dimension (:,:), allocatable :: numx3c_bcna
    integer, dimension (:,:), allocatable :: numy3c_bcna

    real, dimension (:,:), allocatable :: hx_bcna
    real, dimension (:,:), allocatable :: hy_bcna
    real, dimension (:,:), allocatable :: x3cmax_bcna
    real, dimension (:,:), allocatable :: y3cmax_bcna

    real, dimension (:, :, :, :, :), allocatable :: bcna_01
    real, dimension (:, :, :, :, :), allocatable :: bcna_02
    real, dimension (:, :, :, :, :), allocatable :: bcna_03
    real, dimension (:, :, :, :, :), allocatable :: bcna_04
    real, dimension (:, :, :, :, :), allocatable :: bcna_05

    ! XC interactions; 7 implies different derivative types; there are five xc3c arrays - one for each theta
    integer, dimension (:,:), allocatable :: numx3c_xc3c
    integer, dimension (:,:), allocatable :: numy3c_xc3c

    real, dimension (:,:), allocatable :: hx_xc3c
    real, dimension (:,:), allocatable :: hy_xc3c
    real, dimension (:,:), allocatable :: x3cmax_xc3c
    real, dimension (:,:), allocatable :: y3cmax_xc3c

    real, dimension (:, :, :, :, :), allocatable :: xc3c_01
    real, dimension (:, :, :, :, :), allocatable :: xc3c_02
    real, dimension (:, :, :, :, :), allocatable :: xc3c_03
    real, dimension (:, :, :, :, :), allocatable :: xc3c_04
    real, dimension (:, :, :, :, :), allocatable :: xc3c_05
    !xc3c_SN
    ! XC interactions; 7 implies different derivative types; there are five den3 arrays - one for each theta, used only for SNXC method
    integer, dimension (:,:), allocatable :: numx3c_den3
    integer, dimension (:,:), allocatable :: numy3c_den3

    real, dimension (:,:), allocatable :: hx_den3
    real, dimension (:,:), allocatable :: hy_den3
    real, dimension (:,:), allocatable :: x3cmax_den3
    real, dimension (:,:), allocatable :: y3cmax_den3
    real, dimension (:, :, :, :, :), allocatable :: den3_01
    real, dimension (:, :, :, :, :), allocatable :: den3_02
    real, dimension (:, :, :, :, :), allocatable :: den3_03
    real, dimension (:, :, :, :, :), allocatable :: den3_04
    real, dimension (:, :, :, :, :), allocatable :: den3_05
    real, dimension (:, :, :, :, :), allocatable :: den3S_01
    real, dimension (:, :, :, :, :), allocatable :: den3S_02
    real, dimension (:, :, :, :, :), allocatable :: den3S_03
    real, dimension (:, :, :, :, :), allocatable :: den3S_04
    real, dimension (:, :, :, :, :), allocatable :: den3S_05

    !end xc3c_SN

    ! Two center integrals
    ! jel-F2c
    ! integer, dimension (1:21, 0:8) :: ind2c
    integer, dimension (1:23, 0:8) :: ind2c
    ! end jel-F2c
    integer, dimension (:, :, :), allocatable :: numz2c
    real, dimension (:, :, :, :, :), allocatable :: xintegral_2c
    real, dimension (:, :, :, :, :, :), allocatable :: splineint_2c
    real, dimension (:, :, :), allocatable :: z2cmax

    ! One center integrals
    real, dimension (:, :), allocatable :: exc1c_0
    real, dimension (:, :, :, :), allocatable :: exc1c
    real, dimension (:, :, :), allocatable :: xcnu1c
    real, dimension (:, :, :), allocatable :: xcnu1cs
    ! jel-der
    real, dimension (:, :, :), allocatable :: exc1c0
    real, dimension (:, :, :), allocatable :: nuxc1c
    real, dimension (:, :, :, :), allocatable :: dexc1c
    real, dimension (:, :, :), allocatable :: d2exc1c
    real, dimension (:, :, :, :), allocatable :: dnuxc1c
    real, dimension (:, :, :, :, :), allocatable :: d2nuxc1c
    ! end jel-der


    ! AQUI 
    integer, dimension (:), allocatable :: nsu
    integer :: V_intra_dip = 1 ! hacer fdata con vdip_onecenter ? 
    integer, parameter :: nfofx = 207
    integer :: interactions2c_max = 24 !ojo antes estaba en initbasics ufff
    integer, parameter :: numXmax = 31
    integer, parameter :: numYmax = 31
    integer, parameter :: ntheta = 5

contains
    function append_string (filename, extension)
        implicit none
        character (len = 200) append_string
        character (len = 200) filename
        character (len = 200) extension
        integer lenx
        integer leny
        integer len_trim
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
