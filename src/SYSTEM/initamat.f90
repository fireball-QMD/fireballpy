subroutine initamat()
  use M_constants, only: haveDorbitals
  use M_system
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
  amat(:,:,:) = 0.0d0
  if (.not. haveDorbitals) return
  amat(2,1,1) = 0.5d0
  amat(1,2,1) = 0.5d0
  amat(3,2,2) = 0.5d0
  amat(2,3,2) = 0.5d0
  amat(1,1,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
  amat(2,2,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
  amat(3,3,3) = 2.0d0/(2.0d0*sqrt(3.0d0))
  amat(3,1,4) = 0.5d0
  amat(1,3,4) = 0.5d0
  amat(1,1,5) = 0.5d0 
  amat(2,2,5) = - 0.5d0
  return
end

