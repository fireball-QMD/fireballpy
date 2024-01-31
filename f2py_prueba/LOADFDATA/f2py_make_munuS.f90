
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
