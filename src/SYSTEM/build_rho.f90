subroutine build_rho ()
  use iso_c_binding
  use M_system, only: iqout, icluster, igamma, tempfe, blowre, bbnkre, blowim, bbnkim, cape, rho, errno
  implicit none
  rho = 0.0d0
  cape = 0.0d0
  if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
  call denmat ()
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    deallocate (blowre)
    if (igamma .eq. 0) then
      deallocate (blowim)
    end if
  end if
  if (icluster .eq. 0 .and. igamma .eq. 0) deallocate (bbnkim)
  deallocate (bbnkre)

  if (errno .ne. 0) return
  call mixer ()
end subroutine build_rho
