subroutine scf_loop ()
  use M_system
  use M_fdata !, only: nssh !TESTING
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
    write(*,*) 'XXX EBS =',ebs
    write(*,*) 'XXX EBS uiiuee', uiiuee
    write(*,*) 'XXX EBS uxcdcc', uxcdcc
    write(*,*) 'XXX EBS etotxc_1c', etotxc_1c
    write(*,*) 'XXX h_mat =',h_mat(1,1,1,1),h_mat(2,1,1,1)
    write(*,*) 'XXX t_mat =',t_mat(1,1,1,1),t_mat(2,1,1,1)
    write(*,*) 'XXX s_mat =',s_mat(1,1,1,1),s_mat(2,1,1,1)
    write(*,*) 'XXX vna =',vna(1,1,1,1),vna(2,1,1,1)
    write(*,*) 'XXX vxc =',vxc(1,1,1,1),vxc(2,1,1,1)
    write(*,*) 'XXX vxc_1c =',vxc_1c(1,1,1,1),vxc_1c(2,1,1,1)
    write(*,*) 'XXX nuxc1c =',nuxc1c(1,1,1),nuxc1c(2,1,1)
    write(*,*) 'XXX  numorb_max =',numorb_max
    write(*,*) 'XXX  nspecies,nsh_max,nsh_max',nspecies,nsh_max,nsh_max

    stop

    Kscf = Kscf + 1
  end do
end subroutine scf_loop
 
