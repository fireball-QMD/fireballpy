subroutine scf_loop ()
  use M_system
  use M_fdata !, only: nssh !TESTING
  use M_constants
  implicit none
  integer :: iatom,in1,issh !TESTING
  !inicializamos constantes
  scf_achieved = .false.
  delk = 0.0d0
  delk(1,1) = 1.0d0
  delk(2,2) = 1.0d0
  delk(3,3) = 1.0d0
  xlevi = 0.0d0
  xlevi(1,2,3) = 1.0d0
  xlevi(1,3,2) = -1.0d0
  xlevi(3,1,2) = 1.0d0
  xlevi(3,2,1) = -1.0d0
  xlevi(2,3,1) = 1.0d0
  xlevi(2,1,3) = -1.0d0

  Kscf = 1
  do while (Kscf <= max_scf_iterations .and. .not. scf_achieved)
    call assemble_mcweda ()

    call diag_k ()

    call build_rho ()

     write(*,*) 'XXX sigma scf_achieved',sigma, scf_achieved
     do iatom = 1, natoms
      in1 = imass(iatom)
      write (*,'(2x, 10f14.8)') (Qin(issh,iatom), issh = 1, nssh(in1))
      write (*,'(2x, 10f14.8)') (Qout(issh,iatom), issh = 1, nssh(in1))
    end do
    write(*,*)'========== Qout ====== Kscf = ',Kscf
    do iatom = 1, natoms
      in1 = imass(iatom)
      write (*,'(2x, 10f14.8)') (Qout(issh,iatom), issh = 1, nssh(in1))
    end do
    Kscf = Kscf + 1
  end do
end subroutine scf_loop
 
