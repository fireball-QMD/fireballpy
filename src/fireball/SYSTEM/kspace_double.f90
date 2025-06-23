subroutine kspace_double (ikpoint, sks)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: iqout, natoms, degelec, imass, getlssh, getissh, Kscf, blowre, blowim, bbnkre, bbnkim, sm12_complex, eigen_k, norbitals, &
    & norbitals_new, getiatom, neigh_b, neigh_j, neighn, neighPPn, neighPP_b, neighPP_j, vnl, s_mat, h_mat, errno, isgeneig, xl, ratom
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer, intent(in) :: ikpoint
  real(double), dimension(3), intent(in) :: sks
  complex(double), parameter :: a0 = cmplx(0.0d0, 0.0d0, double)
  complex(double), parameter :: a1 = cmplx(1.0d0, 0.0d0, double)
  real(double), parameter :: overtol = 1.0d-4
  integer iatom, imu, inu, in1, in2, ineigh, &
    & jatom, jmu, jnu, mbeta, mineig, lm
  real(double) :: dot
  complex(double) :: phase
  real(double), dimension(3) :: vec
  real(double), dimension (norbitals) :: eigen
  real(double), dimension(3*norbitals-2) :: rwork
  real(double), dimension (:), allocatable :: ww
  complex(double), dimension (:), allocatable :: work
  complex(double), dimension (norbitals, norbitals) :: yyyy, zzzz
  complex(double), dimension (:, :), allocatable :: xxxx
  integer :: lwork, info
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
        vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
        dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
        phase = cmplx(cos(dot), sin(dot), double)
        do inu = 1, num_orb(in2)
          jnu = inu + degelec(jatom)
          do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
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
    sm12_complex(:,:,ikpoint) = zzzz(:,:)
  end if

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
      phase = cmplx(cos(dot), sin(dot), double)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
        end do
      end do
    end do
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
      phase = cmplx(cos(dot), sin(dot), double)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*vnl(imu,inu,ineigh,iatom)
        end do
      end do
    end do
  end do ! do iatom

  if (isgeneig) then
    zzzz(:,:) = sm12_complex(:,:,ikpoint)
    allocate(work(1))
    call zhegv(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, -1, rwork, info)
    lwork = real(work(1))
    deallocate(work)
    allocate(work(lwork))
    call zhegv(1, 'V', 'U', norbitals, yyyy, norbitals, zzzz, norbitals, eigen, work, lwork, rwork, info)
    deallocate(work)
    if (info .eq. 0) then
      norbitals_new = norbitals
      eigen_k(:,ikpoint) = eigen(:)
      bbnkre(:,:,ikpoint) = real(yyyy(:,:))
      bbnkim(:,:,ikpoint) = aimag(yyyy(:,:))
      return
    end if
    if (info .le. norbitals) then
      errno = -info
      return
    end if
    isgeneig = .false.
  end if

  if (Kscf .eq. 1) then
    allocate(work(1))
    call zheev('V', 'U', norbitals, zzzz, norbitals, eigen, work, -1, rwork, info)
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', norbitals, zzzz, norbitals, eigen, work, lwork, rwork, info)
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
      call zgemm('N', 'C', norbitals, norbitals, norbitals, a1, zzzz, norbitals, zzzz, norbitals, a0, xxxx, norbitals)
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
    sm12_complex(:,1:norbitals_new,ikpoint) = xxxx(:,:)
  else
    allocate(xxxx(norbitals, norbitals_new))
    xxxx(:,:) = sm12_complex(:,1:norbitals_new,ikpoint)
  end if
  call zgemm('C', 'N', norbitals_new, norbitals, norbitals, a1, xxxx, norbitals_new, yyyy, norbitals, a0, zzzz(1:norbitals_new,:), norbitals_new)
  call zgemm('N', 'N', norbitals_new, norbitals_new, norbitals, a1, zzzz(1:norbitals_new,:), norbitals_new, xxxx, norbitals, a0, yyyy(1:norbitals_new,1:norbitals_new), norbitals_new)
  allocate(work(1))
  call zheev('V', 'U', norbitals_new, yyyy, norbitals_new, eigen(1:norbitals_new), work, -1, rwork, info)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
  call zheev('V', 'U', norbitals_new, yyyy, norbitals_new, eigen(1:norbitals_new), work, lwork, rwork, info)
  deallocate(work)
  if (info .ne. 0) then
    errno = -info
    return
  end if
  eigen_k(1:norbitals_new,ikpoint) = eigen(1:norbitals_new)
  call zgemm('N', 'N', norbitals, norbitals_new, norbitals_new, a1, xxxx, norbitals, yyyy(1:norbitals_new,1:norbitals_new), norbitals_new, a0, zzzz(:,1:norbitals_new), norbitals)
  bbnkre(:,1:norbitals_new,ikpoint) = real(zzzz(:,1:norbitals_new))
  bbnkim(:,1:norbitals_new,ikpoint) = aimag(zzzz(:,1:norbitals_new))
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    blowre(:,:,ikpoint) = real(yyyy(:,:))
    blowim(:,:,ikpoint) = aimag(yyyy(:,:))
  end if
  deallocate (xxxx)
  if ((Kscf .eq. 1) .and. (iqout .eq. 3)) then
    deallocate (ww)
  end if
end subroutine kspace_double
