subroutine getenergy () 
  use iso_c_binding
  use M_system 
  use M_fdata, only : nssh
  implicit none
  integer(c_long) :: iatom, issh

  call get_ewald (iforce)
  
  call assemble_usr ()
 
  uxcdcc = uxcdcc_ols + etotxc_1c 
  etot = ebs + uiiuee + uxcdcc
  etot = etot + eqmmm
  etotper = etot/natoms
  etotold = etotnew
  etotnew = etotper

end subroutine getenergy

