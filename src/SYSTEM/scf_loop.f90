subroutine scf_loop ()
use M_system
implicit none

do Kscf = 1, max_scf_iterations
  write (*,*) ' Begin scf step = ', Kscf
  call assemble_h ()
  call diag_k ()
  call build_rho ()
  ! mirar scf_achieved en scf_loop_harris 
end do

end subroutine scf_loop
 
