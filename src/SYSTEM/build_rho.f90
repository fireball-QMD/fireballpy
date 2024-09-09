subroutine build_rho ()
  use iso_c_binding
  use M_system
  implicit none
  rho = 0.0d0
  cape = 0.0d0
  if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
  call denmat ()
  call mixer ()
  !if (sigma .lt. sigmatol) scf_achieved = .true.
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
