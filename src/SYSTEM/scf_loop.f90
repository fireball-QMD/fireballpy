subroutine scf_loop ()
use M_system
implicit none

scf_achieved = .false.

do Kscf = 1, max_scf_iterations
  call assemble_mcweda ()
  call diag_k ()
  call build_rho ()
end do

end subroutine scf_loop
 
