subroutine denmat ()
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom
  integer iband
  integer ikpoint
  integer imu, inu
  integer ineigh
  integer in1, in2
  integer iorbital
  integer issh, jssh
  integer jatom
  integer jneigh
  integer mqn
  integer mbeta
  integer mmu
  integer noccupy
  integer nnu
  integer :: info, lwork
  integer, dimension(100) :: work
  real*8 aux1, aux2, aux3
  real*8 deltae
  real*8 dot
  real*8 gutr
  real*8 pcharge
  real*8 ztest
  real*8 checksum
  real*8 Wmu
  real*8, dimension (natoms) :: pqmu
  real*8, dimension (natoms) :: QoutTot
  real*8, dimension (3) :: vec
  complex*16 ai
  complex*16 phase, phasex
  complex*16 step1, step2
  logical read_occupy

  ai = dcmplx(0.0d0,1.0d0)
  rhoPP = 0.0d0
  !AQUI  inquire (file = 'OCCUPATION', exist = read_occupy)

  call fermie ()

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      do ikpoint = 1, nkpoints
        dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) + special_k(3,ikpoint)*vec(3)
        phasex = dcmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
        if (icluster .eq. 0 .and. igamma .eq. 0) then
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*(bbnkre(mmu,iband,ikpoint) - ai*bbnkim(mmu,iband,ikpoint))
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*(bbnkre(nnu,iband,ikpoint) + ai*bbnkim(nnu,iband,ikpoint))
                  gutr = real(step2)
                  rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                  cape(imu,inu,ineigh,iatom) = cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
                end do
              end do
            end if
          end do
        else
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*bbnkre(mmu,iband,ikpoint)
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*bbnkre(nnu,iband,ikpoint)
                  gutr = real(step2)
                  rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                  cape(imu,inu,ineigh,iatom) = cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      do ikpoint = 1, nkpoints
        dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)  + special_k(3,ikpoint)*vec(3)
        phasex = dcmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
        if (icluster .eq. 0 .and. igamma .eq. 0) then
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*(bbnkre(mmu,iband,ikpoint) - ai*bbnkim(mmu,iband,ikpoint))
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*(bbnkre(nnu,iband,ikpoint)  + ai*bbnkim(nnu,iband,ikpoint))
                  gutr = real(step2)
                  rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom) + gutr
                end do
              end do
            end if
          end do
        else
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*bbnkre(mmu,iband,ikpoint)
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*bbnkre(nnu,iband,ikpoint)
                  gutr = real(step2)
                  rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom)  + gutr
                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do

  Qout = 0.0d0
  QLowdin_TOT = 0.0d0

  if (ifixcharge .eq. 1) then
    do iatom = 1, natoms
      in1 = imass(iatom)
      do issh = 1, nssh(in1)
        Qout(issh,iatom) = Qin(issh,iatom)
        QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + Qin(issh,iatom)
      end do
    end do
  else
    if (iqout .eq. 1 .or. iqout .eq. 3) call LOWDIN_CHARGES()

    if (iqout .eq. 2) call MULLIKEN_CHARGES() 

    if (iqout .eq. 4) call MULLIKEN_DIPOLE_CHARGES()

    if (iqout .eq. 7) then
      call MULLIKEN_DIPOLE_CHARGES()
      call Dipole_proyection()
      do iatom = 1, natoms
        in1 = imass(iatom)
        QoutTot(iatom) = 0.0d0
        do imu = 1,nssh(in1)
          QoutTot(iatom) = QoutTot(iatom)+Qout(imu,iatom)
        end do 
      end do 
      do iatom = 1, natoms
        in1 = imass(iatom)
        do imu = 1,nssh(in1)              
          Qout(imu,iatom) = (dq_DP(iatom)/QoutTot(iatom))*Qout(imu,iatom) + Qout(imu,iatom)
        end do 
      end do 
    end if !7
  end if !ifixcharge

  ebs = 0.0d0
  ztest = 0.0d0
  do ikpoint = 1, nkpoints
    do iorbital = 1, norbitals_new
      if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
       ebs = ebs + weight_k(ikpoint)*spin*eigen_k(iorbital,ikpoint)*foccupy(iorbital,ikpoint)
       ztest = ztest + weight_k(ikpoint)*spin*foccupy(iorbital,ikpoint)
      end if
    end do
  end do
  if (abs(ztest - ztot) .gt. 1.0d-02) then
    write (*,*) ' *************** error *************** '
    write (*,*) ' ztest = ', ztest, ' ztot = ', ztot
    write (*,*) ' In denmat.f - ztest .ne. ztot! '
    stop
  end if

end subroutine denmat
