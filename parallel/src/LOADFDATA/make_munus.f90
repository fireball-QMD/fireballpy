subroutine make_munuS ()
  use M_fdata, only: MES_max, nspecies, nssh, index_maxS, muS, mvalueS, nuS
  implicit none
  integer imu, index, in1, in2, issh1, issh2, l1, l2, n1, n2

  muS = 0
  nuS = 0
  mvalueS = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      n1 = 0
      do issh1 = 1 , nssh(in1)
        n1 = n1 + 1
        n2 = 0
        do issh2 = 1, nssh(in2)
          n2 = n2 + 1
          index = index + 1
          muS(index,in1,in2) = n1
          nuS(index,in1,in2) = n2
          mvalueS(index,in1,in2) = 0
        end do
      end do
      index_maxS(in1,in2) = index
    end do
  end do
end subroutine make_munuS
