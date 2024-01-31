! readdata_3c.f90
! Program Description
! ===========================================================================
!       This routine reads the data from the 3-center integral files. When
! read, the information is stored in the array threecint.  This array
! is the field that stores all non-vanishing matrix elements for a general
! 3-center integral.  There are maximal ME3c_max non-vanishing matrix
! elements given on a grid of maximal nfofx data points.  The exact dimensions
! for a given interaction, and a given pair of atoms are numz and num_nonzero.
! ===========================================================================
subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp, maxtype, index, xintegral)
    use M_fdata
    implicit none
    integer, intent (in) :: index
    integer, intent (in) :: iounit
    integer, intent (in) :: isorp
    integer, intent (in) :: maxtype
    integer, intent (in) :: num_nonzero
    integer, intent (in) :: numx, numy 

    real, intent (out), dimension (numXmax, numYmax, ME3c_max, 0:maxtype, nspecies**3) :: xintegral
    !real, intent (inout), dimension (:,:,:,:,:) :: xintegral

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

