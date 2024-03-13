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

    write(*,*)'========== Qin ====== '
    do iatom = 1, natoms
      in1 = imass(iatom)
      write (*,'(2x, 10f14.8)') (Qin(issh,iatom), issh = 1, nssh(in1))
    end do
    write(*,*)'========== Qout ======'
    do iatom = 1, natoms
      in1 = imass(iatom)
      write (*,'(2x, 10f14.8)') (Qout(issh,iatom), issh = 1, nssh(in1))
    end do
    write(*,'(A,F20.6,A,I4,A,F12.10,A,L1)') 'EBS = ',ebs,'; Kscf =',Kscf,'; sigma =',sigma,'; scf_achieved =',scf_achieved
    Kscf = Kscf + 1
  end do
end subroutine scf_loop
 
