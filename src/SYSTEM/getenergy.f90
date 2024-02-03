subroutine getenergy () !MCWEDA only
  use M_system 
  implicit none
  call get_ewald ()
  call assemble_usr (uxcdcc_hf, uiiuee)

  ! to avoid confusion here we add etotxc_1c to double counting term
  ! and set etotxc_1c to zero to do not double it in final print
  uxcdcc = uxcdcc_ols + etotxc_1c 
  etotxc_1c = 0.0d0

  etot = ebs + uiiuee + uxcdcc + etotxc_1c
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
 
