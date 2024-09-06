subroutine scf_loop ()
  use iso_c_binding
  use M_system, only: Kscf,max_scf_iterations,scf_achieved
!  use M_system, only: ebs,sigma
  implicit none
  Kscf = 1
  scf_achieved = .false.
  do Kscf=1,max_scf_iterations+1
    call assemble_mcweda ()
    call diag_k ()
    call build_rho ()
!    write(*,'(3x,A,F15.6,A,I4,A,F12.10)') 'EBS = ',ebs,' ; Kscf =',Kscf,' ; sigma =',sigma
    if(scf_achieved) exit
  end do
end subroutine scf_loop
