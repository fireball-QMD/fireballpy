subroutine calc_me_max (nssh, nsshPP, lssh, lsshPP, &
    & ME2c_max, ME2cDipX_max, ME2cDipY_max, ME2cPP_max, ME3c_max, MES_max)
  use M_fdata, only: nspecies, nsh_max
  implicit none
  integer, dimension(nsh_max), intent(in) :: nssh, nsshPP
  integer, dimension(nsh_max,nspecies), intent(in) :: lssh, lsshPP
  integer, intent(out) :: ME2c_max, ME2cDipX_max, ME2cDipY_max, &
    & ME2cPP_max, ME3c_max, MES_max
  integer index, in1, in2, imu, issh1, issh2, l1, l2, n1, n2

  ! ME2c_max and ME3c_max
  ME2c_max = 0
  ME3c_max = 0
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

  ! ME2cPP_max
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

  ! MES_max
  MES_max = 0
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

  ! ME2cDipX_max and ME2cDipY_max
  ME2cDipX_max = 0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      index = 0
      do issh1 = 1, nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nssh(in2)
          l2 = lssh(issh2,in2)
          if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
            index = index + 1
          end if
          if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
            index = index + 1
          end if
          if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
            index = index + 1
          end if
          if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
            index = index + 1
          end if
          if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
            index = index + 2
          end if
          if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
            index = index + 2
          end if
          if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
            index = index + 1
          end if
          if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
            index = index + 1
          end if
        end do
      end do
      if (index .gt. ME2cDipX_max) ME2cDipX_max = index
    end do
  end do
  ME2cDipY_max = ME2cDipX_max
  if (ME2cDipX_max .gt. ME2c_max) ME2c_max = ME2cDipX_max
end subroutine calc_me_max
