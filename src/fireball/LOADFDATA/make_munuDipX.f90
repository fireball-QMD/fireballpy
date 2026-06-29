subroutine make_munuDipX ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use indices, only: indices_twocenter_set, INDICES_TWOCENTER_DIPX
  use M_fdata, only: ME2cDipX_max, nspecies, nssh, num_orb, lssh, &
    & index_max2cDipX, ME2c_max, muDipX, nuDipX
  implicit none
  integer :: index, in1, in2, issh1, issh2, l1, l2, n1, n2, index_max
  integer, allocatable :: s1(:), s2(:), ll1(:), ll2(:), m1(:), m2(:)

  ME2cDipX_max = 0
  if (allocated(index_max2cDipX)) deallocate(index_max2cDipX)
  allocate (index_max2cDipX (nspecies, nspecies))

  do in1 = 1, nspecies
    num_orb(in1) = 0
    do issh1 = 1 , nssh(in1)
      l1 = lssh(issh1,in1)
      num_orb(in1) = num_orb(in1) + 2*l1 + 1
    end do
  end do

  do in1 = 1, nspecies
    do in2 = 1, nspecies
      call indices_twocenter_set(INDICES_TWOCENTER_DIPX, lssh(1:nssh(in1),in1), &
        &                        lssh(1:nssh(in2),in2), index_max, s1, s2, ll1, ll2, m1, m2)
    end do
    if (index_max > ME2cDipX_max) ME2cDipX_max = index_max
  end do
  if (ME2cDipX_max > ME2c_max) ME2c_max = ME2cDipX_max

  if (allocated(muDipX)) deallocate(muDipX)
  if (allocated(nuDipX)) deallocate(nuDipX)
  allocate (muDipX (ME2cDipX_max, nspecies, nspecies))
  allocate (nuDipX (ME2cDipX_max, nspecies, nspecies))
  muDipX = 0
  nuDipX = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      call indices_twocenter_set(INDICES_TWOCENTER_DIPX, lssh(1:nssh(in1),in1), &
        &                        lssh(1:nssh(in2),in2), index_max, s1, s2, ll1, ll2, m1, m2)
      n1 = 0
      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        n1 = n1 + l1 + 1
        n2 = 0
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          n2 = n2 + l2 + 1
          do index = 1, index_max
            if (issh1 /= s1(index) .or. issh2 /= s2(index) .or. l1 /= ll1(index) .or. l2 /= ll2(index)) cycle
            muDipX(index,in1,in2) = n1 + m1(index)
            nuDipX(index,in1,in2) = n2 + m2(index)
          end do
          n2 = n2 + l2
        end do
        n1 = n1 + l1
      end do
      index_max2cDipX(in1,in2) = index_max
    end do
  end do
end subroutine make_munuDipX
