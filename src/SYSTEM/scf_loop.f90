subroutine scf_loop ()
use M_system
implicit none

do Kscf = 1, max_scf_iterations
  write (*,*) 'Kscf = ', Kscf 
  call assemble_mcweda ()
  write(*,*) 'diag_k'
  call diag_k ()
  write(*,*) 'build_rho'
  call build_rho ()
  write(*,*) 'next Kscf'
  ! mirar scf_achieved en scf_loop_harris 
end do

end subroutine scf_loop
 
