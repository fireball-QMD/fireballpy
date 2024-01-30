! readheader_3c.f90
! Program Description
! ===========================================================================
!       Read the header of the 3-center data files.
!
! ===========================================================================
subroutine readheader_3c (iounit, numx, numy, xmax, ymax)
    use M_fdata
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
