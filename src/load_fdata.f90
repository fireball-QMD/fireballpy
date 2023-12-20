subroutine load_fdata()
  
  use M_fdata
  use M_constants_fireball

  implicit none
  real, dimension (:,:), allocatable :: rcutoff_temp

  open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
  read (12,*) 
  read (12,*) nspecies

  ! TODO read and calculate nsh_max, now is 6 calculate also nssPP_max
  ! TODO lsshPP nsh_max maybe diferent

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

  do ispec = 1, nspecies
    read (12,*)
    read (12,*)
    read (12,102) symbolA(ispec)
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
    read (12,*) (cutoff(issh,ispec), issh = 1, nssh(ispec))
    do issh = 1, nssh(ispec)
     rcutoff(ispec, issh) = rcutoff_temp(issh,ispec)*abohr
    end do
    read (12,103) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
    read (12,103) (napot(issh,ispec), issh = 0, nssh(ispec))
    read (12,*) etotatom(ispec)
    read (12,*)
    ! Jesus borrrrraaaaa..
    if (debug)  then
      write (*,100)
      write (*,301) ispec
      write (*,302) symbolA(ispec)
      write (*,303) nzx(ispec)
      write (*,304) smass(ispec)
      write (*,305) nssh(ispec)
      write (*,306) (lssh(issh,ispec), issh = 1, nssh(ispec))
      write (*,307) nsshPP(ispec)
      write (*,308) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
      write (*,314) rc_PP(ispec)
      write (*,309) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
      write (*,310) (cutoff(issh,ispec), issh = 1, nssh(ispec))
      write (*,311) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
      write (*,312) (napot(issh,ispec), issh = 0, nssh(ispec))
      write (*,313) etotatom(ispec)
      write (*,100)
    endif !debug
  end do !ispec
   
  close(unit = 12) !close info.dat     

  isorpmax = 0
  do in1 = 1, nspecies
    isorpmax = max(isorpmax,nssh(in1))
  end do
  isorpmax_xc = 0
  do in1 = 1, nspecies
    isorpmax_xc = max(isorpmax_xc,nssh(in1))
  end do

  call make_munu (nspecies)
  call make_munuPP (nspecies)
  call make_munuS (nspecies)
  call make_munuDipY (nspecies)
  call make_munuDipX (nspecies)



  ! Procedure progs/READFILES/readdata_mcweda.f90
  ! one-center
  call read_1c (nspecies, itheory, itheory_xc, ispin, ioff2c(7))

  ! two-center
  do interaction = 1, 13
    call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)
  end do
 
  ! Spherical OLSXC exchange-correlation
  do interaction = 15, 23
    call read_2c (interaction, nspecies, itheory, ioff2c(interaction), nzx)
  end do

  interaction = 1   ! bcna
  call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)
  interaction = 3   ! den3 (3c - OLSXC) - average density
  call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)
  interaction = 4   ! den3 (3c - OLSXC) - spherical average density
  call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)

  ! Set up some tables for the 2d interpolator
  call setterp_2d (nspecies, itheory_xc, itheory)

! Deallocate Arrays
! ===========================================================================
  deallocate (rcutoff_temp)
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, a25)
102     format (2x, a2)
103     format (9(2x,a25))
301     format (2x, i2, ' - Information for this species ')
302     format (2x, a2, ' - Element ')
303     format (2x, i3, ' - Nuclear Z ')
304     format (2x, f7.3, ' - Atomic Mass ')
305     format (2x, i2, ' - Number of shells ')
306     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
307     format (2x, i2, ' - Number of shells (Pseudopotential) ')
308     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
309     format (2x, 8(2x,f5.2), ' - Occupation numbers ')
310     format (2x, 8(2x,f5.2), ' - Radial cutoffs ')
311     format (2x, 9(2x,a25), ' - Wavefunction files ')
312     format (2x, 9(2x,a25), ' - (Non)-neutral atom potentials ')
313     format (2x, f12.4, ' - Atomic energy ')
314     format (2x, f12.4, ' - Radial cutoffs (Pseudopotential) ')

end subroutine load_fdata




