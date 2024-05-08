subroutine make_munuPP ()
  use M_fdata, only: nspecies, num_orbPP, lssh, lsshPP, ME2cPP_max, ME2c_max, &
    & index_maxPP, num_orbPP, muPP, nuPP, nssh, nsshPP
  implicit none
  integer imu, index, in1, in2, issh1, issh2, l1, l2, n1, n2

  do in1 = 1, nspecies
    num_orbPP(in1) = 0
    do issh1 = 1 , nsshPP(in1)
      l1 = lsshPP(issh1,in1)
      num_orbPP(in1) = num_orbPP(in1) + 2*l1 + 1
    end do
  end do

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
