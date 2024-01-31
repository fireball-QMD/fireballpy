! make_munuPP.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)) :
!
! num_orbPP (in1) : number of orbitals in pseudoatom-type in1
! muPP (index,in1,in2) : the mu-position for each matrix-element (index)
!                        between pseudoatom-type in1 and pseudoatom-type in2
! nuPP (index,in1,in2) : the nu-position for each matrix-element (index)
!                        between pseudoatom-type in1 and pseudoatom-type in2

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
! (i.e. in MK_INDEX.f)
! ===========================================================================
! Original code written by Jose Ortega.

! Code rewritten by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine make_munuPP ()
    use M_fdata
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
