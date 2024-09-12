subroutine scf_loop (verbose)
  use iso_c_binding
  use M_system, only: max_scf_iterations, scf_achieved, Kscf, etot, sigma, &
    & ebs, eqmmm, uiiuee, uxcdcc_ols, etotxc_1c, sigmatol, iforce, errno
  implicit none
  logical, intent(in) :: verbose
  Kscf = 1
  scf_achieved = .false.
  eqmmm = 0.0d0
  uiiuee = 0.0d0
  uxcdcc_ols = 0.0d0
  etotxc_1c = 0.0d0
  do Kscf=1,max_scf_iterations+1
    call assemble_mcweda ()
    call diag_k ()
    if (errno .ne. 0) return
    call build_rho ()
    if (errno .ne. 0) return

    if (scf_achieved) exit

    if (verbose) then
      ! What is EBS??
      print '(a13,i4,a8,f12.6,a9,e10.3,a3,e10.3)', '-> Iteration ', Kscf, ': ebs = ', ebs, '; debs = ', sigma, ' > ', sigmatol
    end if
  end do
  if (verbose) then
    ! What is EBS??
    print '(a13,i4,a8,f12.6,a9,e10.3,a3,e10.3)', '-> Iteration ', Kscf, ': ebs = ', ebs, '; debs = ', sigma, ' < ', sigmatol
  end if

  call get_ewald (iforce)
  call assemble_usr ()
  etot = ebs + uiiuee + uxcdcc_ols + etotxc_1c + eqmmm
end subroutine scf_loop
