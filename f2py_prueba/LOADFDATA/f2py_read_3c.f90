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
