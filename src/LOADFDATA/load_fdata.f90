subroutine load_fdata()
  use iso_c_binding
  use M_constants, only: abohr
  use M_fdata, only: fdataLocation, infofname, nsh_max, nshPP_max, nspecies, nzx, symbolA, &
    & etotatom, smass, rc_PP, rcutoff, cl_PP, nssh, lssh, nsshPP, lsshPP, Qneutral, wavefxn, &
    & napot, errno2c, errno3c, ind2c, icon3c, splineint_2c, numz2c, z2cmax, interactions2c_max, &
    & isorpmax, isorpmax_xc, ME2c_max, nfofx
  implicit none
  integer(c_long) :: in1, in2, in3, ispec, issh, aux, icount, isorp, &
    & interaction, mintype, maxtype
  real(c_double), dimension (:,:), allocatable :: rcutoff_temp

  ! Find nsh_max and nsh_max_PP
  nsh_max = 0
  nshPP_max = 0
  open (unit = 12, file = trim(fdataLocation) // trim(infofname), status = 'old')
  read (12,*)
  read (12,*) nspecies
  do ispec = 1, nspecies
    do in1 = 1, 5
      read (12,*)
    end do
    read (12,*) aux
    if (aux .gt. nsh_max) nsh_max = aux
    read (12,*)
    read (12,*) aux
    if (aux .gt. nshPP_max) nshPP_max = aux
    do in1 = 1, 8
      read (12,*)
    end do
  end do
  close(12)

  !AQUI pensar, Qinmixer(imix) = Qin(issh,iatom) en miser, (Qinmixer(nsh_max*natoms))
  allocate (nzx (nspecies))
  allocate (symbolA (nspecies))
  allocate (etotatom (nspecies))
  allocate (smass (nspecies))
  allocate (rc_PP (nspecies))
  allocate (rcutoff (nspecies, nsh_max)) 
  allocate (rcutoff_temp (nsh_max, nspecies)) 
  allocate (cl_PP (0:nsh_max - 1, nspecies))
  allocate (nssh (nspecies))
  allocate (lssh (nsh_max, nspecies))
  allocate (nsshPP (nspecies))
  allocate (lsshPP (nsh_max, nspecies))
  allocate (Qneutral (nsh_max, nspecies))
  allocate (wavefxn (nsh_max, nspecies))
  allocate (napot (0:nsh_max, nspecies))
  rcutoff = 0.0d0

  ! Load info.dat
  open (unit = 12, file = trim(fdataLocation) // trim(infofname), status = 'old')
  read (12,*)
  read (12,*)
  do ispec = 1, nspecies
    read (12,*)
    read (12,*)
    read (12,'(2x, a2)') symbolA(ispec)
    read (12,*) nzx(ispec)
    read (12,*) smass(ispec)
    read (12,*) nssh(ispec)
    read (12,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) nsshPP(ispec)
    read (12,*) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
    read (12,*) rc_PP(ispec)
    read (12,*) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
    read (12,*) (rcutoff_temp(issh, ispec), issh = 1, nssh(ispec))
    do issh=1,nssh(ispec)
      rcutoff(ispec,issh) = rcutoff_temp(issh,ispec)*abohr
    end do
    read (12,'(9(2x,a25))') (wavefxn(issh,ispec), issh = 1, nssh(ispec))
    read (12,'(9(2x,a25))') (napot(issh,ispec), issh = 0, nssh(ispec))
    read (12,*) etotatom(ispec)
    read (12,*)
  end do
  close(12)

  ! Compatibility
  isorpmax = nsh_max
  isorpmax_xc = nsh_max

  ! Allocate arrays
  call make_munu ()
  call make_munuPP ()
  call make_munuS ()
  call make_munuDipY ()
  call make_munuDipX ()

  ! Load one centre
  call read_1c ()

  ! Prepare and load two centres
  icount = 0
  ind2c = 0
  errno2c = 0
  interactions2c_max = 27 + 11*nsh_max
  allocate (numz2c (interactions2c_max, nspecies, nspecies))
  allocate (z2cmax (interactions2c_max, nspecies, nspecies))
  allocate (splineint_2c (4, ME2c_max, nfofx, interactions2c_max, nspecies, nspecies))
  numz2c = 0
  z2cmax = 0.0d0
  splineint_2c = 0.0d0
  do interaction = 1, 23
    if ((interaction .ge. 2) .and. (interaction .le. 4)) then
      mintype = 0
      maxtype = nsh_max
    else if ((interaction .ge. 15) .and. (interaction .le. 22)) then
      mintype = 1
      maxtype = nsh_max
    else if ((interaction .ge. 6) .and. (interaction .le. 8)) then
      mintype = 0
      maxtype = 4
    else
      mintype = 0
      maxtype = 0
    end if
    do isorp=mintype,maxtype
      icount = icount + 1
      ind2c(interaction, isorp) = icount
    end do
    if(interaction .eq. 14) cycle
    call read_2c (interaction)
    if (errno2c .ne. 0) return
  end do

  ! Prepare and load three centre
  errno3c = 0
  allocate(icon3c(nspecies, nspecies, nspecies))
  do in1=1,nspecies
    do in2=1,nspecies
      do in3=1,nspecies
        icon3c(in1, in2, in3) = in3 + nspecies*(in2 - 1 + nspecies*(in1 - 1))
      end do
    end do
  end do
  do interaction = 1, 4
    if (interaction .eq. 2) cycle
    call read_3c (interaction)
    if (errno3c .ne. 0) return
  end do
  call setterp_2d ()
  deallocate(rcutoff_temp)
end subroutine load_fdata
