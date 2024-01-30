! make_munuDipY.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)) : (for the case of Y dipole)
!
! num_orb (in1) : number of orbitals in atom-type in1
! mu (index,in1,in2) : the mu-position for each matrix-element (index) between
!                      atom-type in1 and atom-type in2
! nu (index,in1,in2) : the nu-position for each matrix-element (index) between
!                      atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! ***************************************************************************
! JOM: IMPORTANT: in this subroutine, the order of the different matrix
! elements (index) is the SAME as the one given in CREATOR
! (i.e. in mk_indexDipY.f)
! ===========================================================================
subroutine make_munuDipY ()
    use M_fdata
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
