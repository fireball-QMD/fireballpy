subroutine scf_loop (verbose)
  use, intrinsic :: iso_fortran_env, only: dp => real64, stdout => output_unit, stderr => error_unit
  use xc, only: xc_init, xc_end
  use M_fdata, only: exc_names, exc_weights
  use M_system, only: max_scf_iterations, scf_achieved, Kscf, etot, sigma, iqout, &
    & ebs, eqmmm, uiiuee, uxcdcc_ols, etotxc_1c, sigmatol, iforce, errno, Kbest, sigmabest, isgeneig
  implicit none
  logical, intent(in) :: verbose
  Kscf = 1
  scf_achieved = .false.
  isgeneig = .true.
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) isgeneig = .false.
  call xc_init(exc_names, exc_weights)
  do Kscf=1,max_scf_iterations+1
    call assemble_drive()
    if (errno .ne. 0) return
    call diag_k()
    if (errno .ne. 0) return
    call build_rho()
    if (errno .ne. 0) return
    call assemble_usr(0)
    etot = ebs + uiiuee + uxcdcc_ols + etotxc_1c + eqmmm
    if (Kscf .le. max_scf_iterations) call mixer()
    if (scf_achieved .or. (Kscf .gt. max_scf_iterations)) exit
    if (verbose) write (stdout, "('-> i = ',i4,': ETOT = ',f12.6,'; <RMSD_Q> = ',es10.3,' > ',es10.3)") Kscf, etot, sigma, sigmatol
  end do
  if (scf_achieved) then
    if (verbose) write (stdout, "('-> i = ',i4,': ETOT = ',f12.6,'; <RMSD_Q> = ',es10.3,' < ',es10.3)") Kscf, etot, sigma, sigmatol
  else
    write (stderr, "('[WARNING] did not converge, best iteration will be used')")
    if (verbose) write (stdout, "('-> i_best = ',i4,': ETOT = ',f12.6,'; <RMSD_Q> = ',es10.3)") Kbest, etot, sigma
  end if
  call xc_end()
end subroutine scf_loop
