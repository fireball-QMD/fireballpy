subroutine scf_loop ()
  use M_system
  use M_constants
  use M_fdata, only: nssh  
  implicit none
  integer iatom,in1,issh
  Kscf = 1
  scf_achieved = .false.
  do while (Kscf <= max_scf_iterations .and. .not. scf_achieved)
    call assemble_mcweda ()
    call diag_k ()
    call build_rho ()
    Kscf = Kscf + 1
    write(*,'(3x,A,I4,A,F12.10,A,L1)') 'XXX Kscf =',Kscf,'; sigma =',sigma,'; scf_achieved =',scf_achieved
    write (*,'(2x,A, f15.6)') 'XXX EBS = ',ebs
  end do
end subroutine scf_loop
 
