subroutine build_rho ()
  use M_system
  use M_fdata, only: nssh,Qneutral,lssh
  implicit none
  integer i,imu,in1, iatom, matom, ineigh, mbeta, jatom,issh, qmu, ikpoint
  rho = 0.0d0
  cape = 0.0d0
  if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
  call denmat ()
  call mixer ()
  if (sigma .lt. sigmatol) scf_achieved = .true.
  if (sigma .lt. sigmatol) print*,'scf_achieved  Kscf build_rho', scf_achieved,Kscf
  !allocate (blowre_o(norbitals,norbitals,nkpoints))
  if (iqout .ne. 2 .and. icluster .ne. 1) deallocate (blowim)
  if (iqout .ne. 2) deallocate (blowre)
  !deallocate (eigen_k)
  if (icluster .ne. 1) deallocate (bbnkim)
  deallocate (bbnkre)
end subroutine build_rho

