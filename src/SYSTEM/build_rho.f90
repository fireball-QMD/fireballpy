subroutine build_rho ()
  use M_system
  implicit none
  integer imu
  integer ikpoint
  integer i
  rho = 0.0d0
  cape = 0.0d0
  if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
  call denmat (iqout, icluster,tempfe, ebs, bmix, Kscf)
  call mixer (natoms, itheory, ifixcharge, iwrtcharges)
  flag_es = 0 
  if (sigma .lt. sigmatol) then
    scf_achieved = .true.
    flag_es = 1
  endif ! (sigma .lt. sigmatol)
  allocate (bbnkre_o(norbitals,norbitals,nkpoints))
  allocate (blowre_o(norbitals,norbitals,nkpoints))
  if (iqout .ne. 2 .and. icluster .ne. 1) deallocate (blowim)
  if (iqout .ne. 2) deallocate (blowre)
  ! AQUI idynmat iephc ?
  deallocate (eigen_k)
  if (icluster .ne. 1) deallocate (bbnkim)
  deallocate (bbnkre)
  return
end subroutine build_rho

