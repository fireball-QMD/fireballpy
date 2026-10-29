subroutine initamat()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: haveDorbitals
  use M_system, only: amat
  use M_fdata, only: nspecies,lsshPP,nsshPP,lssh,nssh
  implicit none
  integer :: in1
  integer :: issh
  haveDorbitals = .false.
  do in1 = 1, nspecies
   do issh = 1, nssh(in1)
    if (lssh(issh,in1) .eq. 2) haveDorbitals = .true.
   end do
   do issh = 1, nsshPP(in1)
    if (lsshPP(issh,in1) .eq. 2) haveDorbitals = .true.
   end do
  end do
  amat(:,:,:) = 0.0_dp
  if (.not. haveDorbitals) return
  amat(2,1,1) = 0.5_dp
  amat(1,2,1) = 0.5_dp
  amat(3,2,2) = 0.5_dp
  amat(2,3,2) = 0.5_dp
  amat(1,1,3) = - 1.0_dp/(2.0_dp*sqrt(3.0_dp))
  amat(2,2,3) = - 1.0_dp/(2.0_dp*sqrt(3.0_dp))
  amat(3,3,3) = 2.0_dp/(2.0_dp*sqrt(3.0_dp))
  amat(3,1,4) = 0.5_dp
  amat(1,3,4) = 0.5_dp
  amat(1,1,5) = 0.5_dp 
  amat(2,2,5) = - 0.5_dp
  return
end subroutine initamat
