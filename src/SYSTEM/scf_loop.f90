subroutine scf_loop ()
  use M_system
  use M_constants
  use M_fdata, only: nssh 
  implicit none
  Kscf = 1
  do while (Kscf <= max_scf_iterations .and. .not. scf_achieved)
    call assemble_mcweda ()
    call diag_k ()
    call build_rho ()
    write(*,'(A,F20.6,A,I4,A,F12.10,A,L1)') 'EBS = ',ebs,'; Kscf =',Kscf,'; sigma =',sigma,'; scf_achieved =',scf_achieved
    Kscf = Kscf + 1
  end do
end subroutine scf_loop
 
