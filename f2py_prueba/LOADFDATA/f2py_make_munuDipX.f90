

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
