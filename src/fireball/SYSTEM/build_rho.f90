subroutine build_rho ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: iqout, igamma, tempfe, blowre, blowim, cape, rho, errno, Kscf, max_scf_iterations
  implicit none
  rho = 0.0_dp
  cape = 0.0_dp
  if (tempfe .le. 50.0_dp) tempfe = 50.0_dp ! Can't be zero
  call denmat ()
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    deallocate (blowre)
    if (igamma .eq. 0) then
      deallocate (blowim)
    end if
  end if

  if (errno .ne. 0) return
!  if (Kscf .le. max_scf_iterations) call mixer ()
end subroutine build_rho
