subroutine scf_loop ()
  use M_system
  use M_fdata, only: nssh !TESTING
  implicit none
  integer :: iatom,in1,issh !TESTING

  scf_achieved = .false.

  Kscf = 1
  do while (Kscf <= max_scf_iterations .and. .not. scf_achieved)

  call assemble_mcweda ()
  call diag_k ()
  call build_rho ()

  write(*,*)'========== Qout ====== Kscf = ',Kscf
  do iatom = 1, natoms
    in1 = imass(iatom)
    write (*,'(2x, 10f14.8)') (Qout(issh,iatom), issh = 1, nssh(in1))
  end do
  write (*,*)'XXX scf_loop EBS', ebs, Kscf

  Kscf = Kscf + 1
end do

end subroutine scf_loop
 
