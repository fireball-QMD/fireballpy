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
  !allocate (blowre_o(norbitals,norbitals,nkpoints))
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    deallocate (blowre)
    if (igamma .eq. 0) then
      deallocate (blowim)
    end if
  end if
  !deallocate (eigen_k)
  if (icluster .eq. 0 .and. igamma .eq. 0) deallocate (bbnkim)
  deallocate (bbnkre)
end subroutine build_rho

