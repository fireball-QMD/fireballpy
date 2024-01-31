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
