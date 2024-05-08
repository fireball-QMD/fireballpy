subroutine load_fdata()
  use M_fdata, only: fdataLocation, nsh_max, nspecies, &
    & errno2c, errno3c, ind2c, icon3c, xintegral_2c
  implicit none
  integer :: in1, in2, in3, icount, isorp, interaction, mintype, maxtype

  !AQUI pensar, Qinmixer(imix) = Qin(issh,iatom) en miser, (Qinmixer(nsh_max*natoms))

  ! Fill arrays
  call make_munu ()
  call make_munuPP ()
  call make_munuS ()
  call make_munuDipY ()
  call make_munuDipX ()

  ! Load one centre
  call read_1c ()

  ! Prepare and load two centres
  icount = 0
  ind2c = 0
  errno2c = 0
  numz2c = 0
  z2cmax = 0.0d0
  xintegral_2c = 0.0d0
  do interaction = 1, 23
    if ((interaction .ge. 2) .and. (interaction .le. 4)) then
      mintype = 0
      maxtype = nsh_max
    else if ((interaction .ge. 15) .and. (interaction .le. 22)) then
      mintype = 1
      maxtype = nsh_max
    else if ((interaction .ge. 6) .and. (interaction .le. 8)) then
      mintype = 0
      maxtype = 4
    else
      mintype = 0
      maxtype = 0
    end if
    do isorp=mintype,maxtype
      icount = icount + 1
      ind2c(interaction, isorp) = icount
    end do
    if(interaction .eq. 14) cycle
    call read_2c (interaction)
    if (errno2c .ne. 0) return
  end do

  ! Prepare and load three centre
  errno3c = 0
  do in1=1,nspecies
    do in2=1,nspecies
      do in3=1,nspecies
        icon3c(in1, in2, in3) = in3 + nspecies*(in2 - 1 + nspecies*(in1 - 1))
      end do
    end do
  end do
  do interaction = 1, 4
    if (interaction .eq. 2) cycle
    call read_3c (interaction)
    if (errno3c .ne. 0) return
  end do
  call setterp_2d ()
end subroutine load_fdata
