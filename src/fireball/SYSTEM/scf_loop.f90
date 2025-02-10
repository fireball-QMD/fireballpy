subroutine scf_loop (verbose)
  use iso_c_binding
  use M_system, only: max_scf_iterations, scf_achieved, Kscf, etot, sigma, &
    & ebs, eqmmm, uiiuee, uxcdcc_ols, etotxc_1c, sigmatol, iforce, errno, Kbest, sigmabest
  implicit none
  logical, intent(in) :: verbose
  Kscf = 1
  scf_achieved = .false.
  do Kscf=1,max_scf_iterations+1
    call assemble_mcweda ()
    call diag_k ()
    if (errno .ne. 0) return
    call build_rho ()
    if (errno .ne. 0) return
    if (scf_achieved .or. (Kscf .gt. max_scf_iterations)) exit
    if (verbose) print '(a13,i4,a8,f12.6,a9,e10.3,a3,e10.3)', '-> Iteration ', Kscf, ': EBS = ', ebs, '; dEBS = ', sigma, ' > ', sigmatol
  end do
  if (verbose) then
    if (scf_achieved) then
      print '(a13,i4,a8,f12.6,a9,e10.3,a3,e10.3)', '-> Iteration ', Kscf, ': EBS = ', ebs, '; dEBS = ', sigma, ' < ', sigmatol
    else
      print '(a15,i4,a8,f12.6,a9,e10.3)', 'Best iteration ', Kbest, ': EBS = ', ebs, '; dEBS = ', sigmabest
    end if
  end if
  call get_ewald (iforce)
  call assemble_usr ()
  etot = ebs + uiiuee + uxcdcc_ols + etotxc_1c + eqmmm
end subroutine scf_loop
