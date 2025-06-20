subroutine make_munuPP ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nspecies, num_orbPP, lssh, lsshPP, ME2cPP_max, ME2c_max, &
    & index_maxPP, num_orbPP, muPP, nuPP, nssh, nsshPP
  implicit none
  integer :: imu, index, in1, in2, issh1, issh2, l1, l2, n1, n2

  if (allocated(index_maxPP)) deallocate(index_maxPP)
  if (allocated(num_orbPP)) deallocate(num_orbPP)
  allocate (index_maxPP (nspecies, nspecies))
  allocate (num_orbPP (nspecies))
  do in1 = 1, nspecies
    num_orbPP(in1) = 0
    do issh1 = 1 , nsshPP(in1)
      l1 = lsshPP(issh1,in1)
      num_orbPP(in1) = num_orbPP(in1) + 2*l1 + 1
    end do
  end do

  ME2cPP_max = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nsshPP(in2)
          l2 = lsshPP(issh2,in2)
          do imu = -min(l1,l2), min(l1,l2)
            index = index + 1
          end do
        end do
      end do
      if (index .gt. ME2cPP_max) ME2cPP_max = index
    end do
  end do
  if (ME2cPP_max .gt. ME2c_max) ME2c_max = ME2cPP_max

  if (allocated(muPP)) deallocate(muPP)
  if (allocated(nuPP)) deallocate(nuPP)
  allocate (muPP (ME2cPP_max, nspecies, nspecies))
  allocate (nuPP (ME2cPP_max, nspecies, nspecies))
  muPP = 0
  nuPP = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      n1 = 0
      do issh1 = 1 , nssh(in1)
        l1 = lssh(issh1,in1)
        n1 = n1 + l1 + 1
        n2 = 0
        do issh2 = 1, nsshPP(in2)
          l2 = lsshPP(issh2,in2)
          n2 = n2 + l2 + 1
          do imu = -min(l1,l2), min(l1,l2)
            index = index + 1
            muPP(index,in1,in2) = n1 + imu
            nuPP(index,in1,in2) = n2 + imu
          end do
          n2 = n2 + l2
        end do
        n1 = n1 + l1
      end do
      index_maxPP(in1,in2) = index
    end do
  end do
end subroutine make_munuPP
