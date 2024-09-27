subroutine make_munu ()
  use iso_c_binding
  use M_fdata, only: ME2c_max, ME3c_max, num_orb, index_max2c, index_max3c, &
    & nspecies, lssh, nssh, mu, nu, mvalue
  implicit none
  integer(c_long) :: imu, index, in1, in2, issh1, issh2, l1, l2, n1, n2

  ME2c_max = 0
  ME3c_max = 0
  if (allocated(num_orb)) deallocate(num_orb)
  if (allocated(index_max2c)) deallocate(index_max2c)
  if (allocated(index_max3c)) deallocate(index_max3c)
  allocate (num_orb (nspecies))
  allocate (index_max2c (nspecies, nspecies))
  allocate (index_max3c (nspecies, nspecies))
  do in1 = 1, nspecies
    num_orb(in1) = 0
    do issh1 = 1 , nssh(in1)
      l1 = lssh(issh1,in1)
      num_orb(in1) = num_orb(in1) + 2*l1 + 1
    end do
  end do

  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      do issh1 = 1 , nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          do imu = -min(l1,l2), min(l1,l2)
            index = index + 1
          end do
        end do
      end do
      if (index .gt. ME2c_max) ME2c_max = index

      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          if (l1 .eq. 0 .and. l2 .ne. 0) then
            index = index + 1
          end if
          if (l1 .eq. 1) then
            index = index + 1
            if (l2 .ne. 0) then
              index = index + 1
            end if
            if (l2 .eq. 2) then
              index = index + 2
            end if
          end if
          if (l1 .eq. 2) then
            index = index + 1
            if (l2 .ne. 0) then
              index = index + 3
            end if
            if (l2 .eq. 2) then
              index = index + 2
            end if
          end if
        end do
      end do

      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          if (l1 .eq. 2) then
            index = index + 1
          end if
          if (l2 .eq. 2) then
            index = index + 1
          end if
        end do
      end do
      if (index .gt. ME3c_max) ME3c_max = index
    end do
  end do

  if (allocated(mu)) deallocate(mu)
  if (allocated(nu)) deallocate(nu)
  if (allocated(mvalue)) deallocate(mvalue)
  allocate (mu (ME3c_max, nspecies, nspecies))
  allocate (nu (ME3c_max, nspecies, nspecies))
  allocate (mvalue (ME3c_max, nspecies, nspecies))
  mu = 0
  nu = 0
  mvalue = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      n1 = 0
      do issh1 = 1 , nssh(in1)
        l1 = lssh(issh1,in1)
        n1 = n1 + l1 + 1
        n2 = 0
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          n2 = n2 + l2 + 1
          do imu = -min(l1,l2), min(l1,l2)
            index = index + 1
            mu(index,in1,in2) = n1 + imu
            nu(index,in1,in2) = n2 + imu
            mvalue(index,in1,in2) = 0
          end do
        n2 = n2 + l2
        end do
        n1 = n1 + l1
      end do

      index_max2c(in1,in2) = index
      n1 = 0
      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        n1 = n1 + l1 + 1
        n2 = 0
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          n2 = n2 + l2 + 1
          if (l1 .eq. 0 .and. l2 .ne. 0) then
            index = index + 1
            mu(index,in1,in2) = n1
            nu(index,in1,in2) = n2 + 1
            mvalue(index,in1,in2) = 1
          end if
          if (l1 .eq. 1) then
            index = index + 1
            mu(index,in1,in2) = n1 + 1
            nu(index,in1,in2) = n2
            mvalue(index,in1,in2) = 1
            if (l2 .ne. 0) then
              index = index + 1
              mu(index,in1,in2) = n1
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
            end if
            if (l2 .eq. 2) then
              index = index + 1
              mu(index,in1,in2) = n1 + 1
              nu(index,in1,in2) = n2 + 2
              mvalue(index,in1,in2) = 1
              index = index + 1
              mu(index,in1,in2) = n1 - 1
              nu(index,in1,in2) = n2 - 2
              mvalue(index,in1,in2) = 1
            end if
          end if
          if (l1 .eq. 2) then
            index = index + 1
            mu(index,in1,in2) = n1 + 1
            nu(index,in1,in2) = n2
            mvalue(index,in1,in2) = 1
            if (l2 .ne. 0) then
              index = index + 1
              mu(index,in1,in2) = n1
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
              index = index + 1
              mu(index,in1,in2) = n1 - 2
              nu(index,in1,in2) = n2 - 1
              mvalue(index,in1,in2) = 1
              index = index + 1
              mu(index,in1,in2) = n1 + 2
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
            end if
            if (l2 .eq. 2) then
              index = index + 1
              mu(index,in1,in2) = n1 + 1
              nu(index,in1,in2) = n2 + 2
              mvalue(index,in1,in2) = 1
              index = index + 1
              mu(index,in1,in2) = n1 - 1
              nu(index,in1,in2) = n2 - 2
              mvalue(index,in1,in2) = 1
            end if
          end if
          n2 = n2 + l2
        end do
        n1 = n1 + l1
      end do
      n1 = 0
      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        n1 = n1 + l1 + 1
        n2 = 0
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          n2 = n2 + l2 + 1
          if (l1 .eq. 2) then
            index = index + 1
            mu(index,in1,in2) = n1 + 2
            nu(index,in1,in2) = n2
            mvalue(index,in1,in2) = 2
          end if
          if (l2 .eq. 2) then
            index = index + 1
            mu(index,in1,in2) = n1
            nu(index,in1,in2) = n2 + 2
            mvalue(index,in1,in2) = 2
          end if
          n2 = n2 + l2
        end do
        n1 = n1 + l1
      end do
      index_max3c(in1,in2) = index
    end do
  end do
end subroutine make_munu
