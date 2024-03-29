subroutine load_fdata()

    use M_fdata
    use M_constants

    implicit none
    real, dimension (:,:), allocatable :: rcutoff_temp
    integer :: in1
    integer :: ispec
    integer :: issh
    integer :: nsh_max_temp
    integer :: nshPP_max_temp
    integer :: interaction
    integer :: icount, isorp,ideriv

    open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
    read (12,*)
    read (12,*) nspecies

    nsh_max = 0
    nshPP_max = 0
    do ispec = 1, nspecies
        do in1 = 1, 5
            read (12,*)
        end do !in1
            read (12,*) nsh_max_temp
            read (12,*)
            read (12,*) nshPP_max_temp
        do in1 = 1, 8
            read (12,*)
        end do !in1
        if (nsh_max_temp .gt. nsh_max) then
            nsh_max = nsh_max_temp
        end if
        if (nshPP_max_temp .gt. nshPP_max) then
            nshPP_max = nshPP_max_temp
        end if
    end do ! ispec
    close(unit = 12) !close info.dat

    !AQUI pensar, Qinmixer(imix) = Qin(issh,iatom) en miser, (Qinmixer(nsh_max*natoms))
    !nsh_max = 6 !AQUI




    allocate (rcutoff_temp (nsh_max, nspecies))

    allocate (nzx (nspecies))
    allocate (symbolA (nspecies)) 
    allocate (etotatom (nspecies)) 
    allocate (smass (nspecies)) 
    allocate (rc_PP (nspecies)) 
    allocate (rcutoff (nspecies, nsh_max)) 
    rcutoff = 0.0d0
    allocate (cl_PP (0:nsh_max - 1, nspecies))
    allocate (nssh (nspecies))
    allocate (lssh (nsh_max, nspecies))
    allocate (nsshPP (nspecies))
    allocate (lsshPP (nsh_max, nspecies))
    allocate (Qneutral (nsh_max, nspecies))
    allocate (wavefxn (nsh_max, nspecies))
    allocate (napot (0:nsh_max, nspecies))
   
    
    open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
    read (12,*)
    read (12,*)

    do ispec = 1, nspecies
    read (12,*)
    read (12,*)
    read (12,'(2x, a2)') symbolA(ispec)
    read (12,*) nzx(ispec)
    read (12,*) smass(ispec)
    read (12,*) nssh(ispec)
    if (nssh(ispec) .gt. nsh_max) then
        write (*,*) ' nssh(ispec) = ', nssh(ispec),' nsh_max = ', nsh_max
        write (*,*) ' Sorry -- redimension nsh_max in MODULES/dimensions.f90'
        stop
    end if

    read (12,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) nsshPP(ispec)
    read (12,*) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
    read (12,*) rc_PP(ispec)
    read (12,*) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) (rcutoff_temp(issh,ispec), issh = 1, nssh(ispec))
    do issh = 1, nssh(ispec)
    rcutoff(ispec, issh) = rcutoff_temp(issh,ispec)*abohr
    end do
    read (12,'(9(2x,a25))') (wavefxn(issh,ispec), issh = 1, nssh(ispec))
    read (12,'(9(2x,a25))') (napot(issh,ispec), issh = 0, nssh(ispec))
    read (12,*) etotatom(ispec)
    read (12,*)
    end do !ispec

    close(unit = 12) !close info.dat

    isorpmax = 0
    isorpmax_xc = 0
    do in1 = 1, nspecies
    isorpmax = max(isorpmax,nssh(in1))
    isorpmax_xc = max(isorpmax_xc,nssh(in1))
    end do

    icount = 0
    ind2c = 0
    icount = icount + 1
    ind2c(1,0) = icount
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(2,isorp) = icount
    end do
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(3,isorp) = icount
    end do
    do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(4,isorp) = icount
    end do
    icount = icount + 1
    ind2c(5,0) = icount
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(6,ideriv) = icount
    end do
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(7,ideriv) = icount
    end do
    do ideriv = 0, 4
    icount = icount + 1
    ind2c(8,ideriv) = icount
    end do
    icount = icount + 1
    ind2c(9,0) = icount
    icount = icount + 1
    ind2c(10,0) = icount
    icount = icount + 1
    ind2c(11,0) = icount
    icount = icount + 1
    ind2c(12,0) = icount
    icount = icount + 1
    ind2c(13,0) = icount
    icount = icount + 1
    ind2c(14,0) = icount
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(15,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(16,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(17,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(18,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(19,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(20,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(21,isorp) = icount
    end do
    do isorp = 1, isorpmax_xc
    icount = icount + 1
    ind2c(22,isorp) = icount
    end do
    icount = icount + 1
    ind2c(23,0) = icount

    interactions2c_max = icount

    call make_munu ()
    call make_munuPP ()
    call make_munuS ()
    call make_munuDipY ()
    call make_munuDipX ()

    call read_1c ()

    do interaction = 1, 13
    call read_2c (interaction)
    end do

    do interaction = 15, 23
    call read_2c (interaction)
    end do

    interaction = 1   ! bcna
    call read_3c (interaction)
    interaction = 3   ! den3 (3c - OLSXC) - average density
    call read_3c (interaction)
    interaction = 4   ! den3 (3c - OLSXC) - spherical average density
    call read_3c (interaction)

    call setterp_2d ()

    deallocate (rcutoff_temp)

end subroutine load_fdata
