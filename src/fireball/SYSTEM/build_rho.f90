subroutine build_rho ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: iqout, igamma, tempfe, blowre, blowim, cape, rho, errno, Kscf, max_scf_iterations
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

  if (errno .ne. 0) return
  if (Kscf .le. max_scf_iterations) call mixer ()
end subroutine build_rho
