module fdata

    use constants
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

subroutine load_fdata()


    implicit none
    real, dimension (:,:), allocatable :: rcutoff_temp
    integer :: in1
    integer :: ispec
    integer :: issh
    integer :: nsh_max_temp
    integer :: nshPP_max_temp
    integer :: interaction
    integer :: icount, isorp,ideriv

    open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
    read (12,*)
    read (12,*) nspecies

    nsh_max = 0
    nshPP_max = 0
    do ispec = 1, nspecies
    do in1 = 1, 5
    read (12,*)
    end do !in1
    read (12,*) nsh_max_temp
    read (12,*)
    read (12,*) nshPP_max_temp
    do in1 = 1, 8
    read (12,*)
    end do !in1
    if (nsh_max_temp .gt. nsh_max) then
        nsh_max = nsh_max_temp
    end if
    if (nshPP_max_temp .gt. nshPP_max) then
        nshPP_max = nshPP_max_temp
    end if
    end do ! ispec

    ! Not sure if they can be different
    !nsh_max = max(nsh_max, nshPP_max)
    !========================================

    close(unit = 12) !close info.dat

    allocate (rcutoff_temp (nsh_max, nspecies))

    allocate (nzx (nspecies))
    allocate (symbolA (nspecies)) 
    allocate (etotatom (nspecies)) 
    allocate (smass (nspecies)) 
    allocate (rc_PP (nspecies)) 
    allocate (rcutoff (nspecies, nsh_max)) 
    rcutoff = 0.0d0
    allocate (cl_PP (0:nsh_max - 1, nspecies))
    allocate (nssh (nspecies))
    allocate (lssh (nsh_max, nspecies))
    allocate (nsshPP (nspecies))
    allocate (lsshPP (nsh_max, nspecies))
    allocate (Qneutral (nsh_max, nspecies))
    allocate (wavefxn (nsh_max, nspecies))
    allocate (napot (0:nsh_max, nspecies))

    allocate (nsu(nspecies))
    nsu=0 !AQUI

    open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
    read (12,*)
    read (12,*)

    do ispec = 1, nspecies
    read (12,*)
    read (12,*)
    read (12,102) symbolA(ispec)
    read (12,*) nzx(ispec)
    read (12,*) smass(ispec)
    read (12,*) nssh(ispec)
    if (nssh(ispec) .gt. nsh_max) then
        write (*,*) ' nssh(ispec) = ', nssh(ispec),' nsh_max = ', nsh_max
        write (*,*) ' Sorry -- redimension nsh_max in MODULES/dimensions.f90'
        stop
    end if

    read (12,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) nsshPP(ispec)
    read (12,*) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
    read (12,*) rc_PP(ispec)
    read (12,*) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) (rcutoff_temp(issh,ispec), issh = 1, nssh(ispec))
    do issh = 1, nssh(ispec)
    rcutoff(ispec, issh) = rcutoff_temp(issh,ispec)*abohr
    end do
    read (12,103) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
    read (12,103) (napot(issh,ispec), issh = 0, nssh(ispec))
    read (12,*) etotatom(ispec)
    read (12,*)
    ! Jesus borrrrraaaaa..
    if (debug)  then
        write (*,100)
        write (*,301) ispec
        write (*,302) symbolA(ispec)
        write (*,303) nzx(ispec)
        write (*,304) smass(ispec)
        write (*,305) nssh(ispec)
        write (*,306) (lssh(issh,ispec), issh = 1, nssh(ispec))
        write (*,307) nsshPP(ispec)
        write (*,308) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
        write (*,314) rc_PP(ispec)
        write (*,309) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
        write (*,310) (rcutoff_temp(issh,ispec), issh = 1, nssh(ispec))
        write (*,311) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
        write (*,312) (napot(issh,ispec), issh = 0, nssh(ispec))
        write (*,313) etotatom(ispec)
        write (*,100)
    end if !debug
    end do !ispec

    close(unit = 12) !close info.dat

    isorpmax = 0
    isorpmax_xc = 0
    do in1 = 1, nspecies
    isorpmax = max(isorpmax,nssh(in1))
    isorpmax_xc = max(isorpmax_xc,nssh(in1))
    end do



    !ideriv_max = 0
    !if (itheory .eq. 1) ideriv_max = 6

    ! Set up the index field ind2c:
    icount = 0
    ind2c = 0
    icount = icount + 1
    ind2c(1,0) = icount
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(2,isorp) = icount
    end do
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(3,isorp) = icount
    end do
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(4,isorp) = icount
    end do
    icount = icount + 1
    ind2c(5,0) = icount
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(6,ideriv) = icount
    end do
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(7,ideriv) = icount
    end do
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(8,ideriv) = icount
    end do
    icount = icount + 1
    ind2c(9,0) = icount
    icount = icount + 1
    ind2c(10,0) = icount
    icount = icount + 1
    ind2c(11,0) = icount
    icount = icount + 1
    ind2c(12,0) = icount
    icount = icount + 1
    ind2c(13,0) = icount
    icount = icount + 1
    ind2c(14,0) = icount
    !if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4 ) then
    !if (itheory_xc .eq. 4) then 
    !icount = icount + 1
    !ind2c(14,0) = icount
    !end if !end if itheory_xc .eq. 4
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(15,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(16,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(17,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(18,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(19,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(20,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(21,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(22,isorp) = icount
    end do
    icount = icount + 1
    ind2c(23,0) = icount
    !dani.JOM.jel-fr-end


    !end if
    interactions2c_max = icount



    call make_munu ()
    call make_munuPP ()
    call make_munuS ()
    call make_munuDipY ()
    call make_munuDipX ()



    ! Procedure progs/READFILES/readdata_mcweda.f90
    ! one-center

    call read_1c ()

    ! two-center
    do interaction = 1, 13
    call read_2c (interaction)
    end do

    ! interaction = 14 NAC TODO ixczw = 1  

    ! Spherical OLSXC exchange-correlation
    do interaction = 15, 23
    call read_2c (interaction)
    end do

    interaction = 1   ! bcna
    call read_3c (interaction)
    interaction = 3   ! den3 (3c - OLSXC) - average density
    call read_3c (interaction)
    interaction = 4   ! den3 (3c - OLSXC) - spherical average density
    call read_3c (interaction)


    ! Set up some tables for the 2d interpolator
    call setterp_2d ()

    ! Deallocate Arrays
    ! ===========================================================================
    deallocate (rcutoff_temp)

    ! Format Statements
    ! ===========================================================================
    100     format (2x, 70('='))
    101     format (2x, a25)
    102     format (2x, a2)
    103     format (9(2x,a25))
    301     format (2x, i2, ' - Information for this species ')
    302     format (2x, a2, ' - Element ')
    303     format (2x, i3, ' - Nuclear Z ')
    304     format (2x, f7.3, ' - Atomic Mass ')
    305     format (2x, i2, ' - Number of shells ')
    306     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
    307     format (2x, i2, ' - Number of shells (Pseudopotential) ')
    308     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
    309     format (2x, 8(2x,f5.2), ' - Occupation numbers ')
    310     format (2x, 8(2x,f5.2), ' - Radial cutoffs ')
    311     format (2x, 9(2x,a25), ' - Wavefunction files ')
    312     format (2x, 9(2x,a25), ' - (Non)-neutral atom potentials ')
    313     format (2x, f12.4, ' - Atomic energy ')
    314     format (2x, f12.4, ' - Radial cutoffs (Pseudopotential) ')

    return
end subroutine load_fdata

subroutine make_munu ()
    !use dimensions
    !use interactions
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    ! Allocate Arrays
    ! ===========================================================================
    allocate (num_orb (nspecies))
    allocate (index_max2c (nspecies, nspecies))
    allocate (index_max3c (nspecies, nspecies))

    ! Procedure
    ! ===========================================================================
    ! First, calculate the number of orbitals in each atom-type
    do in1 = 1, nspecies
    num_orb(in1) = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    num_orb(in1) = num_orb(in1) + 2*l1 + 1
    end do
    end do

    ! Now we calculate ME3c_max and ME2c_max
    ! ME2c_max is the maximum number of two-center matrix elements: 
    ! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 19
    ! ME3c_max is the maximum number of three-center matrix elements:
    ! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
    ME2c_max = 0
    ME3c_max = 0

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    do imu = -min(l1,l2), min(l1,l2)
    index = index + 1
    end do
    end do
    end do
    if (index .gt. ME2c_max) ME2c_max = index
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)

    if (l1 .eq. 0 .and. l2 .ne. 0) then
        index = index + 1
    end if

    if (l1 .eq. 1) then
        index = index + 1
        if (l2 .ne. 0) then
            index = index + 1
        end if
        if (l2 .eq. 2) then
            index = index + 2
        end if
    end if

    if (l1 .eq. 2) then
        index = index + 1
        if (l2 .ne. 0) then
            index = index + 3
        end if
        if (l2 .eq. 2) then
            index = index + 2
        end if
    end if

    end do
    end do

    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)

    if (l1 .eq. 2) then
        index = index + 1
    end if

    if (l2 .eq. 2) then
        index = index + 1
    end if

    end do
    end do
    if (index .gt. ME3c_max) ME3c_max = index

    end do
    end do

    ! Allocate arrays
    allocate (mu (ME3c_max, nspecies, nspecies))
    allocate (mvalue (ME3c_max, nspecies, nspecies))
    allocate (nu (ME3c_max, nspecies, nspecies))

    ! Now, calculate the mu-nu-map of the matrix elements on the box
    ! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
    ! mu(index) and nu(index).

    ! First, the non-zero two-center (and also three-center) interactions
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    n1 = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    n2 = n2 + l2 + 1
    do imu = -min(l1,l2), min(l1,l2)
    index = index + 1
    mu(index,in1,in2) = n1 + imu
    nu(index,in1,in2) = n2 + imu
    mvalue(index,in1,in2) = 0
    end do
    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do
    index_max2c(in1,in2) = index

    ! For three-center interactions there are some extra non-zero interactions; in
    ! this case, we find out (from symmetry considerations) that non-negative
    ! values of m_i mix with non-negative values of m_j.  Also, negative values of
    ! m_i mix with the negative values of m_j.
    !
    ! Case 1, the interactions with M1 = M2 +- 1      (M=1 case)
    n1 = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    n2 = n2 + l2 + 1

    if (l1 .eq. 0 .and. l2 .ne. 0) then
        index = index + 1
        mu(index,in1,in2) = n1
        nu(index,in1,in2) = n2 + 1
        mvalue(index,in1,in2) = 1
    end if

    if (l1 .eq. 1) then
        index = index + 1
        mu(index,in1,in2) = n1 + 1
        nu(index,in1,in2) = n2
        mvalue(index,in1,in2) = 1

        if (l2 .ne. 0) then
            index = index + 1
            mu(index,in1,in2) = n1
            nu(index,in1,in2) = n2 + 1
            mvalue(index,in1,in2) = 1
        end if

        if (l2 .eq. 2) then
            index = index + 1
            mu(index,in1,in2) = n1 + 1
            nu(index,in1,in2) = n2 + 2
            mvalue(index,in1,in2) = 1

            index = index + 1
            mu(index,in1,in2) = n1 - 1
            nu(index,in1,in2) = n2 - 2
            mvalue(index,in1,in2) = 1
        end if
    end if

    if (l1 .eq. 2) then
        index = index + 1
        mu(index,in1,in2) = n1 + 1
        nu(index,in1,in2) = n2
        mvalue(index,in1,in2) = 1

        if (l2 .ne. 0) then
            index = index + 1
            mu(index,in1,in2) = n1
            nu(index,in1,in2) = n2 + 1
            mvalue(index,in1,in2) = 1

            index = index + 1
            mu(index,in1,in2) = n1 - 2
            nu(index,in1,in2) = n2 - 1
            mvalue(index,in1,in2) = 1

            index = index + 1
            mu(index,in1,in2) = n1 + 2
            nu(index,in1,in2) = n2 + 1
            mvalue(index,in1,in2) = 1
        end if

        if (l2 .eq. 2) then
            index = index + 1
            mu(index,in1,in2) = n1 + 1
            nu(index,in1,in2) = n2 + 2
            mvalue(index,in1,in2) = 1

            index = index + 1
            mu(index,in1,in2) = n1 - 1
            nu(index,in1,in2) = n2 - 2
            mvalue(index,in1,in2) = 1
        end if
    end if
    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do

    ! Case 2, the interactions with M1 = M2 +- 2      (M=2 case)
    n1 = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    n2 = n2 + l2 + 1

    if (l1 .eq. 2) then
        index = index + 1
        mu(index,in1,in2) = n1 + 2
        nu(index,in1,in2) = n2
        mvalue(index,in1,in2) = 2
    end if

    if (l2 .eq. 2) then
        index = index + 1
        mu(index,in1,in2) = n1
        nu(index,in1,in2) = n2 + 2
        mvalue(index,in1,in2) = 2
    end if
    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do
    index_max3c(in1,in2) = index

    ! End loops over the species
    end do
    end do

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine make_munu


subroutine make_munuDipX ()
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    ! Allocate Arrays
    ! ===========================================================================
    allocate (index_max2cDipX (nspecies, nspecies))

    ! Procedure
    ! ===========================================================================
    ! First, calculate the number of orbitals in each atom-type
    do in1 = 1, nspecies
    num_orb(in1) = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    num_orb(in1) = num_orb(in1) + 2*l1 + 1
    end do
    end do


    ! MEDipX_max is the maximum number of two-center matrix elements for X dipole: 

    ME2cDipX_max = 0

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    !         
    if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
        index = index + 1
    end if


    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        index = index + 1
    end if
    !         
    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        index = index + 1
    end if

    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
    end if

    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
    end if

    end do
    end do
    if (index .gt. ME2cDipX_max) ME2cDipX_max = index

    end do
    end do


    ! Allocate arrays
    allocate (muDipX (ME2cDipX_max, nspecies, nspecies))
    allocate (nuDipX (ME2cDipX_max, nspecies, nspecies))

    if (ME2cDipX_max .gt. ME2c_max) ME2c_max = ME2cDipX_max

    ! Now, calculate the mu-nu-map of the matrix elements on the box
    ! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
    ! muDipX(index) and nuDipX(index).

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0 
    n1 = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    n2 = n2 + l2 + 1

    if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1
        nuDipX(index,in1,in2) = n2 + 1
    end if
    !
    if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1 + 1
        nuDipX(index,in1,in2) = n2 
    end if
    !
    if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1
        nuDipX(index,in1,in2) = n2 + 1
    end if
    !
    if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1 + 1
        nuDipX(index,in1,in2) = n2 
    end if


    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1
        nuDipX(index,in1,in2) = n2 + 1

        index = index + 1
        muDipX(index,in1,in2) = n1 + 2
        nuDipX(index,in1,in2) = n2 + 1
    end if
    !
    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1 + 1
        nuDipX(index,in1,in2) = n2 

        index = index + 1
        muDipX(index,in1,in2) = n1 + 1
        nuDipX(index,in1,in2) = n2 + 2
    end if

    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1 - 2
        nuDipX(index,in1,in2) = n2 - 1
    end if

    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipX(index,in1,in2) = n1 - 1
        nuDipX(index,in1,in2) = n2 - 2
    end if

    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do
    index_max2cDipX(in1,in2) = index

    ! end loops over the species
    end do
    end do

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine make_munuDipX

subroutine make_munuDipY ()
    implicit none
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    ! Allocate Arrays
    ! ===========================================================================
    allocate (index_max2cDipY (nspecies, nspecies))

    ! Procedure
    ! ===========================================================================
    ! First, calculate the number of orbitals in each atom-type
    do in1 = 1, nspecies
    num_orb(in1) = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    num_orb(in1) = num_orb(in1) + 2*l1 + 1
    end do
    end do


    ! MEDipY_max is the maximum number of two-center matrix elements for X dipole: 

    ME2cDipY_max = 0

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    !         
    if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
        index = index + 1
    end if
    !         
    if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
        index = index + 1
    end if


    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        index = index + 1
    end if
    !         
    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        index = index + 1
    end if

    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
    end if

    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
    end if

    end do
    end do
    if (index .gt. ME2cDipY_max) ME2cDipY_max = index

    end do
    end do


    ! Allocate arrays
    allocate (muDipY (ME2cDipY_max, nspecies, nspecies))
    allocate (nuDipY (ME2cDipY_max, nspecies, nspecies))

    if (ME2cDipY_max .gt. ME2c_max) ME2c_max = ME2cDipY_max

    ! Now, calculate the mu-nu-map of the matrix elements on the box
    ! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
    ! muDipY(index) and nuDipY(index).

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0 
    n1 = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    l2 = lssh(issh2,in2)
    n2 = n2 + l2 + 1

    if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1
        nuDipY(index,in1,in2) = n2 - 1
    end if
    !
    if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1 - 1
        nuDipY(index,in1,in2) = n2 
    end if
    !
    if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1
        nuDipY(index,in1,in2) = n2 - 1
    end if
    !
    if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1 - 1
        nuDipY(index,in1,in2) = n2 
    end if


    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1
        nuDipY(index,in1,in2) = n2 - 1

        index = index + 1
        muDipY(index,in1,in2) = n1 + 2
        nuDipY(index,in1,in2) = n2 - 1
    end if
    !
    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1 - 1
        nuDipY(index,in1,in2) = n2 

        index = index + 1
        muDipY(index,in1,in2) = n1 - 1
        nuDipY(index,in1,in2) = n2 + 2
    end if

    if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1 - 2
        nuDipY(index,in1,in2) = n2 + 1
    end if

    if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
        index = index + 1
        muDipY(index,in1,in2) = n1 + 1
        nuDipY(index,in1,in2) = n2 - 2
    end if

    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do
    index_max2cDipY(in1,in2) = index

    ! end loops over the species
    end do
    end do

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine make_munuDipY


subroutine make_munuPP ()
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    ! Allocate Arrays
    ! ===========================================================================
    allocate (index_maxPP (nspecies, nspecies))
    allocate (num_orbPP (nspecies))

    ! Procedure
    ! ===========================================================================
    ! First, calculate the number of orbitals in each atom-type
    do in1 = 1, nspecies
    num_orbPP(in1) = 0
    do issh1 = 1 , nsshPP(in1)
    l1 = lsshPP(issh1,in1)
    num_orbPP(in1) = num_orbPP(in1) + 2*l1 + 1
    end do
    end do

    ! Now, calculate ME2cPP_max (JOM); ME2cPP_max is the maximum number of 
    ! two-center matrix elements for PP 
    ME2cPP_max = 0
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    do issh1 = 1, nssh(in1)
    l1 = lssh(issh1,in1)
    do issh2 = 1, nsshPP(in2) 
    l2 = lsshPP(issh2,in2) 
    do imu = -min(l1,l2), min(l1,l2) 
    index = index + 1 
    end do 
    end do 
    end do 
    if (index .gt. ME2cPP_max) ME2cPP_max = index 
    end do 
    end do

    ! Now allocate some arrays (JOM) 
    allocate (muPP (ME2cPP_max, nspecies, nspecies)) 
    allocate (nuPP (ME2cPP_max, nspecies, nspecies)) 

    ! Define max_ME2c, for allocation purposes (JOM)
    if (ME2cPP_max .gt. ME2c_max) ME2c_max = ME2cPP_max

    ! Now, calculate the mu-nu-map of the matrix elements on the box
    ! (num_orb(in1) x num_orbPP(in2)).  For each matrix-element (index) we
    ! calculate muPP(index) and nuPP(index).

    ! First, the non-zero two-center (and also three-center) interactions
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    n1 = 0
    do issh1 = 1 , nssh(in1)
    l1 = lssh(issh1,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do issh2 = 1, nsshPP(in2)
    l2 = lsshPP(issh2,in2)
    n2 = n2 + l2 + 1
    do imu = -min(l1,l2), min(l1,l2)
    index = index + 1
    muPP(index,in1,in2) = n1 + imu
    nuPP(index,in1,in2) = n2 + imu
    end do
    n2 = n2 + l2
    end do
    n1 = n1 + l1
    end do
    index_maxPP(in1,in2) = index

    ! End loops over the species
    end do
    end do

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine make_munuPP

subroutine make_munuS ()
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    ! Allocate Arrays
    ! ===========================================================================
    allocate (index_maxS (nspecies, nspecies))

    ! Procedure
    ! ===========================================================================

    ! Now we calculate MES_max 
    ! MES_max is the maximum number of three-center matrix elements:
    ! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
    MES_max = 0

    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    do issh1 = 1 , nssh(in1)
    do issh2 = 1, nssh(in2)
    index = index + 1
    end do
    end do
    if (index .gt. MES_max) MES_max = index
    end do
    end do

    ! Allocate arrays
    allocate (muS (MES_max, nspecies, nspecies))
    allocate (mvalueS (MES_max, nspecies, nspecies))
    allocate (nuS (MES_max, nspecies, nspecies))

    ! Now, calculate the mu-nu-map of the matrix elements on the box
    ! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
    ! mu(index) and nu(index).

    ! First, the non-zero two-center (and also three-center) interactions
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    index = 0
    n1 = 0
    do issh1 = 1 , nssh(in1)
    n1 = n1 + 1
    n2 = 0
    do issh2 = 1, nssh(in2)
    n2 = n2 + 1
    index = index + 1
    muS(index,in1,in2) = n1
    nuS(index,in1,in2) = n2
    mvalueS(index,in1,in2) = 0
    end do
    end do
    index_maxS(in1,in2) = index

    ! End loops over the species
    end do
    end do

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine make_munuS

subroutine read_1c ()
    implicit none
    integer iline
    integer Nlines_vdip1c_max
    integer trash
    integer in1
    integer ins
    integer issh
    integer isorp
    integer itype
    integer jssh
    integer kssh
    integer kkssh
    integer numsh
    integer, dimension (nsh_max) :: imask
    integer ideriv
    integer iissh, jjssh

    real, dimension (nspecies) :: idshell

    logical skip_it

    character (len=70) message
    character (len = 200) extension
    character (len = 200) filename
    character (len = 200) root
    character(2) :: auxz    

    !*******************************************************************
    !          M c W E D A   E X C H A N G E - C O R R E L A T I O N
    !******************************************************************* 
    !if (itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then 

    allocate(exc1c0 (nspecies,nsh_max,nsh_max))
    allocate(nuxc1c (nspecies,nsh_max,nsh_max))
    allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
    allocate(d2exc1c (nspecies,nsh_max,nsh_max))
    allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
    allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))

    !==========================================================
    !  READ FILE             nxc1c_dqi.XX.dat
    !==========================================================
    do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/xc1c_dqi.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
    read (36,100) message
    end do
    read (36,*) itype, numsh
    do issh = 1, numsh
    read (36,*) (exc1c0(in1,issh,jssh),jssh=1,numsh)
    end do
    read (36,*)
    do issh = 1, numsh
    read (36,*) (nuxc1c(in1,issh,jssh),jssh=1,numsh)
    end do
    ! exc1c0(in1) = exc1c0(in1)
    do issh = 1, numsh
    do jssh = 1, numsh
    nuxc1c(in1,issh,jssh) = nuxc1c(in1,issh,jssh)
    exc1c0(in1,issh,jssh) = exc1c0(in1,issh,jssh)
    end do
    end do
    close(36)
    end do !in1 .. nspecies

    !===========================================================
    !  READ FILE             nuxc1crho.XX.dat
    !===========================================================

    do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/nuxc1crho.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
    read (36,100) message
    end do
    read (36,*) itype, numsh, kkssh

    do kssh = 1, nssh(in1)
    do issh = 1, numsh
    read (36,*) (dnuxc1c(in1,issh,jssh,kssh),jssh=1,numsh)
    end do
    do issh = 1, numsh
    do jssh = 1, numsh
    dnuxc1c(in1,issh,jssh,kssh) = dnuxc1c(in1,issh,jssh,kssh)
    end do
    end do
    end do ! do kssh
    close(36)
    end do !in1 .. nspecies

    !=========================================================
    !  READ FILE             exc1crho.XX.dat
    !=========================================================

    do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/exc1crho.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
    read (36,100) message
    end do
    read (36,*) itype, numsh, kkssh
    do kssh = 1, nssh(in1)
    do issh = 1, numsh
    read (36,*) (dexc1c(in1,issh,jssh,kssh),jssh=1,numsh)
    end do
    do issh = 1, numsh
    do jssh = 1, numsh
    dexc1c(in1,issh,jssh,kssh) = dexc1c(in1,issh,jssh,kssh)
    end do
    end do
    end do ! do kssh
    close(36)
    end do !in1 .. nspecies

    !if (itheory_xc .eq. 4) then
    !end if !end if itheory_xc .eq. 4
    !end if ! if(itheory_xc.eq.2 .or. itheory_xc .eq. 4) 

    !+++++++++++++++++++++++++++++++NEW JUNE 2019+++++++++++++++++++++++++++
    !.........................Vip 1c...........................................
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (V_intra_dip .eq. 1) then
        allocate(Nlines_vdip1c(nspecies))
        Nlines_vdip1c_max=0
        root = trim(fdataLocation)//'/vdip_onecenter'
        do in1 = 1,nspecies
        write (extension,'(''_'',i2.2)') in1
        filename = append_string (root,extension)
        open (unit = 36, file = filename, status = 'unknown')
        read(36,*) Nlines_vdip1c(in1)
        if (Nlines_vdip1c(in1) .gt. Nlines_vdip1c_max) then
            Nlines_vdip1c_max=Nlines_vdip1c(in1)
        end if
        close(36)
        end do !end do in1

        allocate(muR(Nlines_vdip1c_max,nspecies))
        allocate(nuR(Nlines_vdip1c_max,nspecies))
        allocate(alphaR(Nlines_vdip1c_max,nspecies))
        allocate(betaR(Nlines_vdip1c_max,nspecies))
        allocate(IR(Nlines_vdip1c_max,nspecies))
        muR    = 0.0d0
        nuR    = 0.0d0
        alphaR = 0.0d0
        betaR  = 0.0d0
        IR     = 0.0d0

        do in1 = 1,nspecies
        write (extension,'(''_'',i2.2)') in1
        filename = append_string (root,extension)
        open (unit = 36, file = filename, status = 'unknown')
        read(36,*) trash
        do iline = 1,Nlines_vdip1c(in1)
        read(36,*) muR(iline,in1), nuR(iline,in1), alphaR(iline,in1), betaR(iline,in1), IR(iline,in1)
        end do !end do iline = 1,Nlines_vdip1c
        close(36)
        end do !end do in1 = 1,nspecies
    end if ! if (V_intra_dip .eq. 1)

    100     format (a70)
    200     format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
    return
end subroutine read_1c

subroutine read_2c (interaction)
    implicit none

    integer, intent (in) :: interaction
    integer in1
    integer in2
    integer initype
    integer iounit
    integer isorp
    integer issh
    integer isub2c
    integer itype
    integer maxtype
    integer npseudo
    integer num_nonzero
    integer numz
    integer nzx1
    integer nzx2

    real rc1
    real rc2
    real zmax
    real zmin

    ! We include up to 5 non-local (L values) of the pseudopotential. Usually,
    ! you will have 2 (L=0, 1 and sometimes 2 (D)).
    real, dimension (nsh_max) :: cl_pseudo

    !character (len = 200) append_string
    character (len = 200) extension
    character (len = 200) filename
    character (len = 200) root
    character (len = 200) root_isorp

    !external append_string


    ! Procedure
    ! ===========================================================================
    ! Set up isub2c. 
    isub2c = 0
    if (interaction .eq. 6) isub2c = 4
    if (interaction .eq. 7) isub2c = 4
    if (interaction .eq. 8) isub2c = 4

    ! Allocate and initialize arrays
    if(interaction .eq. 1) then
        write(*,*)ME2c_max, nfofx, interactions2c_max,nspecies, nspecies
        allocate (xintegral_2c (ME2c_max, nfofx, interactions2c_max,nspecies, nspecies))
        allocate (z2cmax (interactions2c_max, nspecies, nspecies))
        allocate (numz2c (interactions2c_max, nspecies, nspecies))
        xintegral_2c = 0.0d0
    end if
    iounit = 71

    ! Here are the roots for the file names for all interactions
    if (interaction .eq. 1)  root = trim(fdataLocation)//'/overlap'
    if (interaction .eq. 2)  root = trim(fdataLocation)//'/vna_ontopl'
    if (interaction .eq. 3)  root = trim(fdataLocation)//'/vna_ontopr'
    if (interaction .eq. 4)  root = trim(fdataLocation)//'/vna_atom  '
    if (interaction .eq. 5)  root = trim(fdataLocation)//'/vnl       '
    if (interaction .eq. 6)  root = trim(fdataLocation)//'/xc_ontop'
    if (interaction .eq. 7)  root = trim(fdataLocation)//'/xc_atom   '
    if (interaction .eq. 8)  root = trim(fdataLocation)//'/xc_corr   '
    if (interaction .eq. 9)  root = trim(fdataLocation)//'/dipole_z  '
    if (interaction .eq. 10) root = trim(fdataLocation)//'/dipole_y  '
    if (interaction .eq. 11) root = trim(fdataLocation)//'/dipole_x  '
    if (interaction .eq. 12) root = trim(fdataLocation)//'/coulomb   '
    if (interaction .eq. 13) root = trim(fdataLocation)//'/kinetic   '
    if (interaction .eq. 14) root = trim(fdataLocation)//'/nuxc      '
    if (interaction .eq. 15) root = trim(fdataLocation)//'/den_ontopl'
    if (interaction .eq. 16) root = trim(fdataLocation)//'/den_ontopr'
    if (interaction .eq. 17) root = trim(fdataLocation)//'/den_atom  ' 
    if (interaction .eq. 18) root = trim(fdataLocation)//'/dnuxc_ol  '
    if (interaction .eq. 19) root = trim(fdataLocation)//'/dnuxc_or  '
    if (interaction .eq. 20) root = trim(fdataLocation)//'/denS_ontopl'
    if (interaction .eq. 21) root = trim(fdataLocation)//'/denS_ontopr'
    if (interaction .eq. 22) root = trim(fdataLocation)//'/denS_atom  ' 
    if (interaction .eq. 23) root = trim(fdataLocation)//'/overlapS  ' 

    ! Now generate the file name of the file to be opened.  Loop over all cases of
    ! the interaction (e.g. different charges of the xc-stuff)




    ! Loop over atoms in1
    do in1 = 1, nspecies

    ! Loop over atoms in2
    do in2 = 1, nspecies

    ! Initialize initype and maxtype
    initype = 0
    if(interaction .ge. 15 ) initype = 1
    if(interaction .eq. 23 ) initype = 0
    maxtype = 0
    if (interaction .eq. 2) isub2c = nssh(in1)
    if (interaction .eq. 3) isub2c = nssh(in2)
    if (interaction .eq. 4) isub2c = nssh(in2)
    if (interaction .eq. 15) isub2c = nssh(in1)
    if (interaction .eq. 16) isub2c = nssh(in2)
    if (interaction .eq. 17) isub2c = nssh(in2)
    if (interaction .eq. 18) isub2c = nssh(in1)
    if (interaction .eq. 19) isub2c = nssh(in2)
    if (interaction .eq. 20) isub2c = nssh(in1)
    if (interaction .eq. 21) isub2c = nssh(in2)
    if (interaction .eq. 22) isub2c = nssh(in2)

    !if (itheory .eq. 1) maxtype = isub2c
    maxtype = isub2c
    ! Harris case for average density
    if (interaction .ge. 15 .and. interaction .le. 22) maxtype = isub2c 
    do isorp = initype, maxtype

    ! Append the number of subtypes on root if there is more than one file for a
    ! given pair of atoms and a given interaction.  The result will be root_isorp
    root_isorp = root
    if (isub2c .ge. 1) then    
        write (extension,'(''_'',i2.2)') isorp
        root_isorp = append_string (root,extension)
    end if

    ! Append the nuclear charges at root, the result will be root_isorp.nz1.nz2.dat
    nzx1 = nzx(in1)
    nzx2 = nzx(in2)
    write (extension,'(''.'',i2.2,''.'',i2.2)') nzx1, nzx2
    filename = append_string (root_isorp, extension)
    write (extension,'(''.dat'')')
    filename = append_string (filename, extension)
    !           if (isorp .eq. initype) write (*,'('' Opening data file: '',a100)') filename
    open (unit = iounit, file = filename, status = 'old')
    call readheader_2c (interaction, iounit, nsh_max, numz, rc1, rc2, &
        &                         zmin, zmax, npseudo, cl_pseudo)
    if (numz .gt. nfofx) then
        write (*,*) ' numz = ', numz, ' in read_2c.f90'
        write (*,*) ' nfofx = ',nfofx
        write (*,*) ' Fix this parameter and recompile! '
        stop
    end if

    ! For the Non-local pseudopotential ONLY.
    if (interaction .eq. 5) then
        do issh = 1, nsshPP(in2) 
        cl_PP(issh,in2) = cl_pseudo(issh)
        end do
    end if

    ! Here are the data file characteristics: number of points and the grid range
    itype = ind2c(interaction,isorp)
    z2cmax(itype,in1,in2) = zmax
    numz2c(itype,in1,in2) = numz

    ! The array index_max2c(in1,in2) is the number of non-vanishing matrix elements
    ! for a general 2-center integral that are stored in the field
    ! twocint(index_max2c(in1,in2),numz,j2x,in1,in2) at numz different bond charge
    ! distances.
    num_nonzero = index_max2c(in1,in2)

    ! For pseudopotential.
    if (interaction .eq. 5) num_nonzero = index_maxPP(in1,in2)

    ! For the vna_atom and xc_atoms cases, the number of interactions is equal
    ! to the number of shells for in1,in1 pair not in1, in2 pair.  This is because
    ! the wavefunctions are both located only on in1!
    if (interaction .eq. 4 .or. interaction .eq. 7)                   &
        &      num_nonzero = index_max2c(in1,in1)
    ! JIMM
    if (interaction .eq. 10)  num_nonzero = index_max2cDipY(in1,in2)
    if (interaction .eq. 11)  num_nonzero = index_max2cDipX(in1,in2)
    !
    if (interaction .eq. 17)  num_nonzero = index_max2c(in1,in1)
    if (interaction .eq. 22)  num_nonzero = index_maxS(in1,in1)

    ! Special case for coulomb (short-range) part.  It's special because it is not
    ! a matrix element.  This in interaction = 12
    if (interaction .eq. 12 .or. interaction .eq. 14)                 &
        &      num_nonzero = nssh(in1)*nssh(in2)
    ! Special case for spehrical density approximation.  It's special because 
    ! it has different size of a matrix element.  
    ! This in interaction = 20,21,22
    if (interaction .eq. 20 .or. interaction .eq. 21)                 &
        &       num_nonzero = index_maxS(in1,in2)
    if (interaction .eq. 23) num_nonzero = index_maxS(in1,in2)

    call readdata_2c (interaction, iounit, num_nonzero, numz, zmax,itype, in1, in2)

    close (unit = iounit)
    end do
    end do ! nspecies
    end do ! nspecies

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ========================================================================


    return
end subroutine read_2c
subroutine read_3c (interaction)
    implicit none
    integer, intent (in) :: interaction
    integer iounit
    integer in1, in2, in3
    integer index
    integer isorp
    integer itheta
    integer maxtype, mintype
    integer numx, numy
    integer nz1, nz2, nz3

    real xmax
    real ymax

    !character (len=200) append_string
    character (len=200) extension
    character (len=200) filename
    character (len=200) root
    character (len=200) root1
    character (len=200) root2

    !external append_string

    ! Allocate Arrays
    ! ===========================================================================

    ! Procedure
    ! ===========================================================================
    ! Allocate and initialize integral arrays to zero
    if (interaction .eq. 1) then
        maxtype = isorpmax
        allocate (bcna_01(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (bcna_02(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (bcna_03(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (bcna_04(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (bcna_05(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (numx3c_bcna(0:maxtype, nspecies**3))
        allocate (numy3c_bcna(0:maxtype, nspecies**3))
        allocate (hx_bcna(0:maxtype, nspecies**3))
        allocate (hy_bcna(0:maxtype, nspecies**3))
        allocate (x3cmax_bcna(0:maxtype, nspecies**3))
        allocate (y3cmax_bcna(0:maxtype, nspecies**3))

        bcna_01 = 0.0d0
        bcna_02 = 0.0d0
        bcna_03 = 0.0d0
        bcna_04 = 0.0d0
        bcna_05 = 0.0d0
        numx3c_bcna = 0
        numy3c_bcna = 0
        x3cmax_bcna = 0.0d0
        y3cmax_bcna = 0.0d0
    else if (interaction .eq. 2) then
        maxtype = ideriv_max
        allocate (xc3c_01(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (xc3c_02(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (xc3c_03(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (xc3c_04(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (xc3c_05(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (numx3c_xc3c(0:maxtype, nspecies**3))
        allocate (numy3c_xc3c(0:maxtype, nspecies**3))
        allocate (hx_xc3c(0:maxtype, nspecies**3))
        allocate (hy_xc3c(0:maxtype, nspecies**3))
        allocate (x3cmax_xc3c(0:maxtype, nspecies**3))
        allocate (y3cmax_xc3c(0:maxtype, nspecies**3))

        xc3c_01 = 0.0d0
        xc3c_02 = 0.0d0
        xc3c_03 = 0.0d0
        xc3c_04 = 0.0d0
        xc3c_05 = 0.0d0
        numx3c_xc3c = 0
        numy3c_xc3c = 0
        x3cmax_xc3c = 0.0d0
        y3cmax_xc3c = 0.0d0

    else if (interaction .eq. 3) then

        maxtype = isorpmax_xc
        allocate (den3_01(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3_02(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3_03(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3_04(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3_05(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (numx3c_den3(0:maxtype, nspecies**3))
        allocate (numy3c_den3(0:maxtype, nspecies**3))
        allocate (hx_den3(0:maxtype, nspecies**3))
        allocate (hy_den3(0:maxtype, nspecies**3))
        allocate (x3cmax_den3(0:maxtype, nspecies**3))
        allocate (y3cmax_den3(0:maxtype, nspecies**3))

        den3_01 = 0.0d0
        den3_02 = 0.0d0
        den3_03 = 0.0d0
        den3_04 = 0.0d0
        den3_05 = 0.0d0
        numx3c_den3 = 0
        numy3c_den3 = 0
        x3cmax_den3 = 0.0d0
        y3cmax_den3 = 0.0d0

    else if (interaction .eq. 4 ) then 
        ! sphere approx.
        maxtype = isorpmax_xc
        allocate (den3S_01(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3S_02(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3S_03(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3S_04(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        allocate (den3S_05(numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3))
        den3S_01 = 0.0d0
        den3S_02 = 0.0d0
        den3S_03 = 0.0d0
        den3S_04 = 0.0d0
        den3S_05 = 0.0d0

    end if

    ! Set unit number for opening data file.
    iounit = 71

    ! Set up the index field icon3c:
    if (.not. allocated(icon3c)) then
        allocate (icon3c (nspecies, nspecies, nspecies))
        do in1 = 1, nspecies
        do in2 = 1, nspecies
        do in3 = 1, nspecies
        icon3c(in1,in2,in3) = in3 + (in2-1)*nspecies + (in1-1)*nspecies**2
        end do
        end do
        end do
    end if

    ! Generate the file name of the file to be opened
    ! Here are the roots of the file names for all interactions
    if (interaction .eq. 1) root1 = trim(fdataLocation)//'/bcna'
    if (interaction .eq. 2) root1 = trim(fdataLocation)//'/xc3c'
    if (interaction .eq. 3) root1 = trim(fdataLocation)//'/den3'
    if (interaction .eq. 4) root1 = trim(fdataLocation)//'/deS3'

    ! Loop over all species in1, in2, in3
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    do in3 = 1, nspecies
    index = icon3c(in1,in2,in3)

    ! Loop over all cases of the interaction (e.g. different charges of the 
    ! xc-stuff, or the different contributions for the neutral atom potential stuff)
    maxtype = 0
    mintype = 0
    !if (itheory .eq. 1) then
    if (interaction .eq. 1) maxtype = nssh(in3)
    if (interaction .eq. 2) maxtype = ideriv_max
    !end if
    ! Average density: we need loop over all orbitals on ialp atoms           
    if (interaction .eq. 3 .or. interaction .eq. 4) then 
        maxtype = nssh(in3)
        mintype = 1
    end if

    do isorp = mintype, maxtype

    ! Loop over all expansion coefficients
    do itheta = 1, ntheta

    ! Append the number of expansion coefficients which is the number of files for 
    ! each interaction subtype.  The result will be root_itheta.
    write(extension,'(''_'',i2.2)') itheta
    root2 = append_string(root1,extension)
    root = root2

    ! Append the number of subtypes on root if there is more than one interaction 
    ! subtype.  The result will be root_itheta_isorp.
    write(extension,'(''_'',i2.2)') isorp
    root = append_string(root2,extension)

    ! Append the atom numbers at root, the result will be 
    ! root_itheta_isorp.nz1.nz2.nz3.dat
    nz1 = nzx(in1)
    nz2 = nzx(in2)
    nz3 = nzx(in3)
    write(extension,'(''.'',i2.2,''.'',i2.2,''.'',i2.2)') nz1, nz2, nz3
    filename = append_string(root,extension)
    write(extension,'(''.dat'')')
    filename = append_string(filename,extension)

    !             if (itheta .eq. 1 .and. isorp .eq. 0)                           &
    !     &        write(*,'('' Opening data file: '',a100)')filename
    !             if (itheta .eq. 1 .and. isorp .eq. 1 .and. interaction .eq. 3)  &
    !     &        write(*,'('' Opening data file: '',a100)') filename
    !             if (itheta .eq. 1 .and. isorp .eq. 1 .and. interaction .eq. 4)  &
    !     &        write(*,'('' Opening data file: '',a100)') filename
    open (unit = iounit, file = filename, status = 'old')


    ! ymin, ymax, numy: bond charge distances, grid
    ! xmin, xmax, numx: neutral atom distances, grid
    call readheader_3c (iounit, numx, numy, xmax, ymax)

    ! Here are the data file characteristics: number of points and the grid range
    if (itheta .eq. 1) then
        if (interaction .eq. 1) then
            x3cmax_bcna(isorp,index) = xmax
            y3cmax_bcna(isorp,index) = ymax
            numx3c_bcna(isorp,index) = numx
            numy3c_bcna(isorp,index) = numy
        else if (interaction .eq. 2) then
            x3cmax_xc3c(isorp,index) = xmax
            y3cmax_xc3c(isorp,index) = ymax
            numx3c_xc3c(isorp,index) = numx
            numy3c_xc3c(isorp,index) = numy
        else if (interaction .eq. 3 .or. interaction .eq. 4) then
            x3cmax_den3(isorp,index) = xmax
            y3cmax_den3(isorp,index) = ymax
            numx3c_den3(isorp,index) = numx
            numy3c_den3(isorp,index) = numy
        end if
    end if

    ! The variable inter_max3c(in1,in2) is the number of non-vanishing matrix
    ! elements for a general 3-center integral that are stored in the fields
    ! bcna(inter_max3c(in1,in2),numx,numy,0:numorb_max,icon3c(in1,in2,in3)) and
    ! xc3c(inter_max3c(in1,in2),numx,numy,0:iderivmax,icon3c(in1,in2,in3)) 
    ! We have to specify the atoms, the exact interaction type, and the number of 
    ! the expansion coefficient.
    if (interaction .eq. 1) then
        if (itheta .eq. 1) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax,  index,  bcna_01)
        if (itheta .eq. 2) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax,  index,  bcna_02)
        if (itheta .eq. 3) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax,  index,  bcna_03)
        if (itheta .eq. 4) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax,  index,  bcna_04)
        if (itheta .eq. 5) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax,  index,  bcna_05)
    else if (interaction .eq. 2 ) then
        if (itheta .eq. 1) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,ideriv_max,index,  xc3c_01)
        if (itheta .eq. 2) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,ideriv_max,index,  xc3c_02)
        if (itheta .eq. 3) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,ideriv_max,index,  xc3c_03)
        if (itheta .eq. 4) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,ideriv_max,index,  xc3c_04)
        if (itheta .eq. 5) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,ideriv_max,index,  xc3c_05)
    else if (interaction .eq. 3) then

        if (itheta .eq. 1) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax_xc, index, den3_01)
        if (itheta .eq. 2) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax_xc, index, den3_02)
        if (itheta .eq. 3) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax_xc, index, den3_03)
        if (itheta .eq. 4) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax_xc, index, den3_04)
        if (itheta .eq. 5) call readdata_3c (iounit, numx, numy, index_max3c(in1,in2),isorp,isorpmax_xc, index, den3_05)
    else if (interaction .eq. 4) then 
        if (itheta .eq. 1) call readdata_3c (iounit, numx, numy, index_maxS(in1,in2), isorp,isorpmax_xc,index, den3S_01)
        if (itheta .eq. 2) call readdata_3c (iounit, numx, numy, index_maxS(in1,in2), isorp,isorpmax_xc,index, den3S_02)
        if (itheta .eq. 3) call readdata_3c (iounit, numx, numy, index_maxS(in1,in2), isorp,isorpmax_xc,index, den3S_03)
        if (itheta .eq. 4) call readdata_3c (iounit, numx, numy, index_maxS(in1,in2), isorp,isorpmax_xc,index, den3S_04)
        if (itheta .eq. 5) call readdata_3c (iounit, numx, numy, index_maxS(in1,in2), isorp,isorpmax_xc,index, den3S_05)
    end if
    close (iounit)
    end do   !  interaction subtype
    end do     !   \
    end do       !     -----  loops over in1,in2,in3
    end do         !   /
    end do           !   expansion coefficients loop

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine read_3c
subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
    implicit none
    integer, intent (in) :: in1, in2
    integer, intent (in) :: interaction
    integer, intent (in) :: iounit
    integer, intent (in) :: itype
    integer, intent (in) :: num_nonzero
    integer, intent (in) :: numz
    real, intent (in) :: zmax
    integer ipoint
    integer integral

    real, dimension (ME2c_max, nfofx) :: gstore

    ! Allocate Arrays
    ! ===========================================================================

    ! Procedure
    ! ===========================================================================


    if (interaction .ne. 8) then
        do ipoint = 1, numz
        read (iounit,*) (gstore(integral,ipoint), integral = 1, num_nonzero)
        end do
        do ipoint = 1, numz
        do integral = 1, num_nonzero
        xintegral_2c(integral,ipoint,itype,in1,in2) =  gstore(integral,ipoint)
        end do
        end do
    else
        do ipoint = 1, numz
        read (iounit,*) gstore(1,ipoint)
        xintegral_2c(1,ipoint,itype,in1,in2) = gstore(1,ipoint)
        end do
    end if
    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end subroutine readdata_2c
subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
    implicit none
    integer, intent (in) :: index
    integer, intent (in) :: iounit
    integer, intent (in) :: isorp
    integer, intent (in) :: maxtype
    integer, intent (in) :: num_nonzero
    integer, intent (in) :: numx, numy 

    !real, intent (out), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3) :: xintegral
    real, intent (inout), dimension (:,:,:,:,:) :: xintegral

    ! Local Parameters and Data Declaration
    ! ===========================================================================

    ! Local Variable Declaration and Description
    ! ===========================================================================
    integer ipoint
    integer integral
    integer jpoint

    real, dimension (ME3c_max, numXmax, numYmax) :: gstore
    real, allocatable, save, dimension(:,:) :: binomial
    integer maxmax
    real, external :: factorial

    ! Procedure
    ! ===========================================================================
    do jpoint = 1, numy
    do ipoint = 1, numx
    read (iounit,*) (gstore(integral,ipoint,jpoint), integral = 1, num_nonzero)
    end do
    end do

    do jpoint = 1, numy
    do ipoint = 1, numx
    do integral = 1, num_nonzero
    xintegral(ipoint,jpoint,integral,isorp,index) = gstore(integral,ipoint,jpoint)  
    end do
    end do
    end do

    return
end subroutine readdata_3c






subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1,   &
        &                            rc2, zmin, zmax, npseudo, cl_pseudo)
    implicit none

    ! Argument Declaration and Description
    ! ===========================================================================
    ! Input
    integer, intent (in) :: interaction
    integer, intent (in) :: iounit
    integer, intent (in) :: nsh_max

    ! Output
    integer, intent (out) :: npseudo
    integer, intent (out) :: numz

    real, intent (out) :: rc1, rc2
    real, intent (out) :: zmin, zmax

    real, intent (out), dimension (nsh_max) :: cl_pseudo

    ! Local Parameters and Data Declaration
    ! ===========================================================================

    ! Local Variable Declaration and Description
    ! ===========================================================================
    integer iline
    integer issh
    integer nucz1, nucz2

    character (len = 70) message

    ! Allocate Arrays
    ! ===========================================================================

    ! Procedure
    ! ===========================================================================
    do iline = 1, 9
    read (iounit,100) message
    end do

    read (iounit,*) nucz1, rc1
    read (iounit,*) nucz2, rc2

    ! The pseudopotential has an extra 2 lines.
    if (interaction .eq. 5) then
        read (iounit,*) npseudo
        read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
    end if

    read (iounit,*) zmax, numz
    zmin = 0.0d0

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================
    100     format (a70)

    return
end subroutine readheader_2c
subroutine readheader_3c (iounit, numx, numy, xmax, ymax)
    implicit none
    integer, intent (in) :: iounit
    integer, intent (out) :: numx, numy
    real, intent (out) :: xmax, ymax
    integer iline
    integer nucZ1, nucZ2, nucZ3, nr, ntheta_in, nphi2
    real rc1a, rc2a, rc3a

    character (len = 70) message

    ! Allocate Arrays
    ! ===========================================================================

    ! Procedure
    ! ===========================================================================
    do iline = 1, 10
    read (iounit,100) message
    end do

    read (iounit,*) nphi2, nr, ntheta_in

    read (iounit,*) ymax, numy
    read (iounit,*) xmax, numx

    read (iounit,100) message

    read (iounit,*) nucZ1, rc1a
    read (iounit,*) nucZ2, rc2a
    read (iounit,*) nucZ3, rc3a

    read (iounit,100) message

    if (numx .gt. numXmax .or. numy .gt. numYmax) then
        write (*,*) ' Courseness too fine in 3c data files. '
        write (*,*) ' numx = ', numx, ' numXmax = ', numXmax
        write (*,*) ' numy = ', numy, ' numYmax = ', numYmax
        write (*,*) ' Change numXmax and numYmax in MODULES/dimensions.f90! '
        stop
    end if

    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================
    100     format (a70)

    return
end subroutine readheader_3c
subroutine setterp_2d ()
    implicit none
    integer in1
    integer in2
    integer in3
    integer index
    integer isorp

    real xmin
    real ymin

    ! Procedure
    ! ===========================================================================
    write (*,*) ' Running setterp_2d. Set up two-dimensional interpolator. '

    ! Table some trivial things that otherwise get recomputed at every possible
    ! opportunity.

    ! For now, we assume ALWAYS that xmin = 0.0d0
    xmin = 0.0d0
    ymin = 0.0d0
    do in1 = 1, nspecies
    do in2 = 1, nspecies
    do in3 = 1, nspecies
    index = icon3c(in1,in2,in3)
    do isorp = 0, isorpmax
    hx_bcna(isorp,index) = (x3cmax_bcna(isorp,index) - xmin)         &
        &                            /(numx3c_bcna(isorp,index) - 1)
    hy_bcna(isorp,index) = (y3cmax_bcna(isorp,index) - ymin)         &
        &                            /(numy3c_bcna(isorp,index) - 1)
    end do
    !if (itheory .ne. 3) then 
    !if(itheory_xc .eq. 0) then
    !do isorp = 0, ideriv_max
    !hx_xc3c(isorp,index) = (x3cmax_xc3c(isorp,index) - xmin)  &
    !&                            /real(numx3c_xc3c(isorp,index) - 1)
    !hy_xc3c(isorp,index) = (y3cmax_xc3c(isorp,index) - ymin)  &
    !&                            /real(numy3c_xc3c(isorp,index) - 1)
    !end do
    !else if(itheory_xc .ne. 3) then
    do isorp = 1, isorpmax_xc
    hx_den3(isorp,index) = (x3cmax_den3(isorp,index) - xmin)   &
        &                            /real(numx3c_den3(isorp,index) - 1)
    hy_den3(isorp,index) = (y3cmax_den3(isorp,index) - ymin)   &
        &                            /real(numy3c_den3(isorp,index) - 1)
    end do
    !end if
    !endif

    end do   ! end do in3
    end do  ! end do in2 
    end do ! end do in1

    ! Format Statements
    ! ===========================================================================

    return
end subroutine setterp_2d
end module fdata
