subroutine getenergy () 
  use iso_c_binding
  use M_system 
  implicit none

  call get_ewald (iforce)
  call assemble_usr ()
 
  uxcdcc = uxcdcc_ols + etotxc_1c 
  etot = ebs + uiiuee + uxcdcc
  etot = etot + eqmmm
  etotper = etot/natoms
  etotold = etotnew
  etotnew = etotper

end subroutine getenergy
