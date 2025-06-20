subroutine kspace_gamma ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: iqout, natoms, degelec, imass, getlssh, getissh, Kscf, blowre, bbnkre, sm12_real, eigen_k, norbitals, &
    & norbitals_new, getiatom, neigh_b, neigh_j, neighn, neighPPn, neighPP_b, neighPP_j, vnl, s_mat, h_mat, errno, isgeneig
  use M_fdata, only: num_orb, Qneutral
  implicit none
  real(double), parameter :: a0 = 0.0d0
  real(double), parameter :: a1 = 1.0d0
  real(double), parameter :: overtol = 1.0d-4
  integer ikpoint, iatom, imu, inu, in1, in2, ineigh, &
    & jatom, jmu, jnu, mbeta, mineig, lm
  real(double), dimension (norbitals) :: eigen
  real(double), dimension (:), allocatable :: ww, work
  real(double), dimension (:, :), allocatable :: xxxx, yyyy, zzzz
  integer :: liwork, lwork, info
  integer, dimension (:), allocatable :: iwork
  ikpoint = 1
  allocate(yyyy(norbitals,norbitals))
  allocate(zzzz(norbitals,norbitals))
  yyyy = a0

  if ((Kscf .eq. 1) .and. (iqout .eq. 3)) then
    allocate(ww(norbitals))
    do inu = 1, norbitals
      imu = getissh(inu)
      iatom = getiatom(inu)
      lm = getlssh(inu)
      in1 = imass(iatom)
      if(Qneutral(getissh(inu),imass(getiatom(inu))).lt.0.01)then
        ww(inu) = 1.0d0
      else
        ww(inu) = 10.0d0
      end if
    end do
  end if

  if (Kscf .eq. 1) then
    zzzz = a0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
        mbeta = neigh_b(ineigh,iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        do inu = 1, num_orb(in2)
          jnu = inu + degelec(jatom)
          do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + s_mat(imu,inu,ineigh,iatom)
          end do
        end do
      end do
    end do ! do iatom
    if (iqout .eq. 3) then
      do inu = 1, norbitals
        do imu = 1, norbitals
          zzzz(inu,imu) = zzzz(inu,imu)*ww(inu)*ww(imu)
        end do
      end do
    end if
    sm12_real(:,:) = zzzz(:,:)
  end if

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + h_mat(imu,inu,ineigh,iatom)
        end do
      end do
    end do
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + vnl(imu,inu,ineigh,iatom)
        end do
      end do
    end do
  end do ! do iatom

  if (isgeneig) then
    zzzz(:,:) = sm12_real(:,:)
    allocate(xxxx(norbitals,norbitals))
    xxxx(:,:) = yyyy(:,:)
    allocate(work(1))
    call dsygv(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, -1, info)
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call dsygv(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, lwork, info)
    deallocate(work)
    !allocate(work(1), iwork(1))
    !call dsygvd(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, -1, iwork, -1, info)
    !lwork = work(1)
    !liwork = iwork(1)
    !deallocate(work, iwork)
    !allocate(work(lwork))
    !allocate(iwork(liwork))
    !call dsygv(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, lwork, iwork, liwork, info)
    !deallocate(work, iwork)
    if (info .eq. 0) then
      norbitals_new = norbitals
      eigen_k(:,ikpoint) = eigen(:)
      bbnkre(:,:,ikpoint) = real(yyyy(:,:), double)
      deallocate (xxxx)
      deallocate (yyyy)
      deallocate (zzzz)
      return
    end if
    if (info .le. norbitals) then
      write (*,*) 'me muero'
      open (unit=360, file='ham.dat', status='unknown')
      do inu = 1, norbitals
        do imu = 1, norbitals
      write(360,*) xxxx(imu,inu), zzzz(imu,inu)
      end do
      end do
      close(360)
      errno = -info
      return
    end if
    isgeneig = .false.
  end if

  if (Kscf .eq. 1) then
    allocate(work(1))
    call dsyev('V', 'U', norbitals, zzzz, norbitals, eigen, work, -1, info)
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call dsyev('V', 'U', norbitals, zzzz, norbitals, eigen, work, lwork, info)
    deallocate(work)
    if (info .ne. 0) then
      errno = -info
      return
    end if
    mineig = 0
    do imu = 1, norbitals
      if (eigen(imu) .lt. overtol) mineig = imu
    end do
    mineig = mineig + 1
    norbitals_new = norbitals + 1 - mineig
    if ((norbitals_new .ne. norbitals) .and. ((iqout .eq. 1) .or. (iqout .eq. 3))) then
      errno = 13
      return
    end if
    allocate(xxxx(norbitals, norbitals_new))
    if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
      do imu = 1, norbitals
        zzzz(:,imu) = zzzz(:,imu)*eigen(imu)**(-0.25d0)
      end do
      call dgemm('N', 'C', norbitals, norbitals, norbitals, a1, zzzz, norbitals, zzzz, norbitals, a0, xxxx, norbitals)
      if (iqout .eq. 3) then
        do imu=1, norbitals
          xxxx(imu,:)=xxxx(imu,:)*ww(imu)
        end do
      end if
    else
      do imu = mineig, norbitals
        xxxx(:,imu-mineig+1) = zzzz(:,imu)*eigen(imu)**(-0.5d0)
      end do
    end if
    sm12_real(:,1:norbitals_new) = xxxx(:,:)
  else
    allocate(xxxx(norbitals, norbitals_new))
    xxxx(:,:) = sm12_real(:,1:norbitals_new)
  end if
  call dgemm('C', 'N', norbitals_new, norbitals, norbitals, a1, xxxx, norbitals_new, yyyy, norbitals, a0, zzzz(1:norbitals_new,:), norbitals_new)
  call dgemm('N', 'N', norbitals_new, norbitals_new, norbitals, a1, zzzz(1:norbitals_new,:), norbitals_new, xxxx, norbitals, a0, yyyy(1:norbitals_new,1:norbitals_new), norbitals_new)
  allocate(work(1))
  call dsyev('V', 'U', norbitals_new, yyyy, norbitals_new, eigen(1:norbitals_new), work, -1, info)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
  call dsyev('V', 'U', norbitals_new, yyyy, norbitals_new, eigen(1:norbitals_new), work, lwork, info)
  deallocate(work)
  if (info .ne. 0) then
    errno = -info
    return
  end if
  eigen_k(1:norbitals_new,ikpoint) = eigen(1:norbitals_new)
  call dgemm('N', 'N', norbitals, norbitals_new, norbitals_new, a1, xxxx, norbitals, yyyy(1:norbitals_new,1:norbitals_new), norbitals_new, a0, zzzz(:,1:norbitals_new), norbitals)
  bbnkre(:,1:norbitals_new,ikpoint) = real(zzzz(:,1:norbitals_new), double)
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) blowre(:,:,ikpoint) = real(yyyy(:,:), double)
  deallocate (xxxx)
  deallocate (yyyy)
  deallocate (zzzz)
  if ((Kscf .eq. 1) .and. (iqout .eq. 3)) then
    deallocate (ww)
  end if
end subroutine kspace_gamma
