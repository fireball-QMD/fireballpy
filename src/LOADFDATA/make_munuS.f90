subroutine make_munuS ()
  use iso_c_binding
  use M_fdata, only: MES_max, nspecies, nssh, index_maxS, muS, mvalueS, nuS
  implicit none
  integer(c_long) :: index, in1, in2, issh1, issh2, n1, n2

  MES_max = 0
  if (allocated(index_maxS)) deallocate(index_maxS)
  allocate (index_maxS (nspecies, nspecies))
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      do issh1 = 1 , nssh(in1)
        do issh2 = 1, nssh(in2)
          index = index + 1
        end do
      end do
      if (index .gt. MES_max) MES_max = index
    end do
  end do

  if (allocated(muS)) deallocate(muS)
  if (allocated(nuS)) deallocate(nuS)
  if (allocated(mvalueS)) deallocate(mvalueS)
  allocate (muS (MES_max, nspecies, nspecies))
  allocate (nuS (MES_max, nspecies, nspecies))
  allocate (mvalueS (MES_max, nspecies, nspecies))
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
