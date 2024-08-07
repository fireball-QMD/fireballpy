subroutine mixer ()
  use M_system
  use M_fdata, only: nssh, nsh_max
  implicit none
  integer iatom
  integer in1
  integer issh
  integer imix
  real*8 dqrms
  real*8 dqmax
  real*8 renorm
  real*8 zcheck
  real*8 zouttot
  dqrms = 0
  dqmax = -99
  Qinmixer=0.0d0
  Qoutmixer=0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      dqmax = amax1(abs(Qin(issh,iatom) - Qout(issh,iatom)),dqmax)
      dqrms = dqrms + (Qin(issh,iatom) - Qout(issh,iatom))**2
    end do
  end do
  dqrms = sqrt(dqrms)/(2*natoms)
  imix = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      imix = imix + 1
      Qinmixer(imix) = Qin(issh,iatom)
      Qoutmixer(imix) = Qout(issh,iatom)
    end do
  end do

  call anderson (Qoutmixer, Qinmixer, imix)

  if (Kscf .gt. 1) then
    if (sigma .lt. sigmaold) then
      sigmaold = sigma
    end if
  else
    sigmaold = sigma
  end if

  !if (sigma .lt. sigmatol) then
  !  scf_achieved = .true.
  !  deallocate(Fv, Xv, delF, delX, r2_sav)
  !end if
  if (.not. scf_achieved) then
    imix = 0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do issh = 1, nssh(in1)
        imix = imix + 1
        Qin(issh,iatom) = Qinmixer(imix)
      end do
    end do
  end if
  zouttot = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      zouttot = zouttot + Qin(issh,iatom)
    end do
  end do
  renorm = (zouttot - ztot)/nssh_tot
  zcheck = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Qin(issh,iatom) = Qin(issh,iatom) - renorm
      zcheck = zcheck + Qin(issh,iatom)
    end do
  end do
end subroutine mixer

