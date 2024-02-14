subroutine read_1c ()
    use M_fdata
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

    allocate(exc1c0 (nspecies,nsh_max,nsh_max))
    allocate(nuxc1c (nspecies,nsh_max,nsh_max))
    allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
    allocate(d2exc1c (nspecies,nsh_max,nsh_max))
    allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
    allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))

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
        do issh = 1, numsh
            do jssh = 1, numsh
                nuxc1c(in1,issh,jssh) = nuxc1c(in1,issh,jssh)
                exc1c0(in1,issh,jssh) = exc1c0(in1,issh,jssh)
            end do
        end do
        close(36)
    end do !in1 .. nspecies

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
end subroutine read_1c

