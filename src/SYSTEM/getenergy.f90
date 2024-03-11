subroutine getenergy () 
  use M_system 
  use M_fdata, only : nssh
  implicit none
  integer :: iatom, issh

  call get_ewald (iforce)
  
  call assemble_usr ()
 
  uxcdcc = uxcdcc_ols + etotxc_1c 
  etot = ebs + uiiuee + uxcdcc
  etot = etot + eqmmm
  etotper = etot/natoms

  write (*,*) ' ---------- T H E  T O T A L  E N E R G Y ----------- '
  write (*,*) '  '
  write (*,502) ebs
  write (*,503) uiiuee
  write (*,504) etotxc_1c
  write (*,505) uxcdcc
  write (*,507) etot
  write (*,508) etotper
  write (*,509) atomic_energy
  write (*,510) etot - atomic_energy
  write (*,512) efermi
  write (*,*) '  '
  write (*,511) (etot - atomic_energy)/natoms
  write (*,*) ' ----------------------------------------------------- '

  etotold = etotnew
  etotnew = etotper

  write(*,*)'========== Qout ======'
  do iatom = 1, natoms
    write (*,'(2x, 10f14.8)') (Qout(issh,iatom), issh = 1, nssh(imass(iatom)))
  end do


  ! Format Statements
  ! =======================================================
100     format (2x, 70('='))
500     format (2x, ' Time step = ', i6, ' SCF step = ', i3)
501     format (2x, ' Time step = ', i6)
502     format (2x, '           ebs = ', f15.6)
503     format (2x, '     uii - uee = ', f15.6)
504     format (2x, '     etotxc_1c = ', f15.6)
505     format (2x, '        uxcdcc = ', f15.6)
507     format (2x, '          ETOT = ', f15.6)
508     format (2x, '     Etot/atom = ', f15.6)
509     format (2x, ' Atomic Energy = ', f15.6)
510     format (2x, '     CohesiveE = ', f15.6)
511     format (2x, ' Cohesive Energy per atom  = ', f15.6)
512     format (2x, '   Fermi Level = ', f15.6)

end subroutine getenergy

