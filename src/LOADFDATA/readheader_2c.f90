subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1,   &
        &                            rc2, zmin, zmax, npseudo, cl_pseudo)
    implicit none
    integer, intent (in) :: interaction
    integer, intent (in) :: iounit
    integer, intent (in) :: nsh_max
    integer, intent (out) :: npseudo
    integer, intent (out) :: numz
    real, intent (out) :: rc1, rc2
    real, intent (out) :: zmin, zmax
    real, intent (out), dimension (nsh_max) :: cl_pseudo
    integer iline
    integer issh
    integer nucz1, nucz2
    character (len = 70) message

    do iline = 1, 9
        read (iounit,100) message
    end do
    read (iounit,*) nucz1, rc1
    read (iounit,*) nucz2, rc2

    if (interaction .eq. 5) then
        read (iounit,*) npseudo
        read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
    end if

    read (iounit,*) zmax, numz
    zmin = 0.0d0

    100     format (a70)

    return
end subroutine readheader_2c
