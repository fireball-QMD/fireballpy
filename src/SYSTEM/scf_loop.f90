subroutine scf_loop ()
  use M_system,only : Kscf,max_scf_iterations,scf_achieved
  use M_fdata, only : nssh  
  implicit none
  integer iatom,in1,issh
  Kscf = 1
  scf_achieved = .false.
  do while (Kscf <= max_scf_iterations .and. .not. scf_achieved)
    call assemble_mcweda ()
    call diag_k ()
    call build_rho ()
    Kscf = Kscf + 1
  end do
end subroutine scf_loop
 
