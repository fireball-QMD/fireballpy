!AQUI ufff
subroutine kspace ( ikpoint, sks)
  use M_system
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer, intent (in) :: ikpoint
  real, intent (in), dimension (3) :: sks
  real*8, parameter :: overtol = 1.0d-4
  integer iatom
  integer imu
  integer info
  integer inu
  integer in1, in2
  integer ineigh
  integer ishort
  integer jatom
  integer jmu
  integer jnu
  integer mbeta
  integer mineig
  integer lm
  integer issh
  real sqlami
  real*8, dimension (norbitals) :: eigen
  real*8, dimension (norbitals) :: slam
  real*8 a0
  real*8 a1
  real*8 magnitude
  real*8, dimension (:, :), allocatable :: xxxx
  real*8, dimension (:, :), allocatable :: yyyy
  real*8, dimension (:, :), allocatable :: zzzz
  real*8, dimension (:, :), allocatable, save :: sm12_save
  real*8, dimension (:, :), allocatable :: ssss
  real*8, dimension (:), allocatable :: ww
  real*8, allocatable, dimension (:) :: work
  integer, allocatable, dimension (:) :: iwork
  integer lwork, liwork
  magnitude = sqrt(sks(1)**2 + sks(2)**2 + sks(3)**3)
  if (magnitude .gt. 1.0d-3) then
    write (*,*) ' gamma point, the complex version of LAPACK is needed.'
    write (*,*) ' We must stop here! '
    stop
  end if
  a0 = 0.0d0
  a1 = 1.0d0
  ishort = 1
  allocate (xxxx(norbitals,norbitals))
  allocate (yyyy(norbitals,norbitals))
  allocate (zzzz(norbitals,norbitals))
  if (iqout .eq. 3) then
   allocate (ssss(norbitals,norbitals))
   allocate (ww(norbitals))
  endif
  liwork = 15*norbitals
  allocate (iwork(liwork))
  lwork = 1
  allocate(work(lwork))
  if (.not. allocated(sm12_save)) then
    allocate (sm12_save(norbitals,norbitals))
  end if
  zzzz = a0
  yyyy = a0
  xxxx = a0
  if (iqout .eq. 3) then
    do inu = 1, norbitals
      imu = getissh(inu)
      iatom = getiatom(inu)
      lm = getlssh(inu)
      in1 = imass(iatom)
      if(Qneutral(getissh(inu),imass(getiatom(inu))).lt.0.01)then
        ww(inu) = 1.0d0
      else
        ww(inu) = 10.0d0
      endif
    enddo
  endif

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
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + h_mat(imu,inu,ineigh,iatom)
        end do ! do imu
      end do ! do inu
    end do ! do ineigh
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      do inu = 1, num_orb(in2)
        jnu = inu + degelec(jatom)
        do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          yyyy(jmu,jnu) = yyyy(jmu,jnu) + vnl(imu,inu,ineigh,iatom)
        end do ! do imu
      end do ! do inu
    end do ! do inegh
  end do ! do iatom
  if (iqout .eq. 3) then
    do inu = 1, norbitals
      do imu = 1, norbitals
        ssss(inu,imu) = zzzz(inu,imu)
      end do
    end do
    do inu = 1, norbitals
      do imu = 1, norbitals
        zzzz(inu,imu) = zzzz(inu,imu)*ww(inu)*ww(imu)
      end do
    end do
  endif  ! end if (iqout .eq. 3)
  if (Kscf .eq. 1 .or. iqout .eq. 3) then
    call dsyev ('V', 'U', norbitals, zzzz, norbitals, slam, work, -1, info)
    lwork = work(1)
    deallocate (work)
    allocate(work(lwork))
    call dsyev ('V', 'U', norbitals, zzzz, norbitals, slam, work,lwork, info )
    if (info .ne. 0) call diag_error (info, 0)
    mineig = 0
    do imu = 1, norbitals
      if (slam(imu) .lt. overtol) mineig = imu
    end do
    mineig = mineig + 1
    norbitals_new = norbitals + 1 - mineig
    if (norbitals_new .ne. norbitals) then
      write (*,*) ' Linear dependence encountered in basis set. '
      if(ishort .eq. 1) then    ! Don't print out again if done above
        write (*,*) '      The overlap eigenvalues: '
      else          ! They asked for extra printout
        write(*,*) ' Eigenvectors that correspond to eigenvalues'
        do imu = 1, mineig - 1
          write(*,*) ' eigenvector',imu
          do jmu = 1, norbitals
            write(*,*) jmu,' ',zzzz(jmu,imu)
          end do
        end do
      end if
      do imu = mineig, norbitals
        jmu = imu - mineig + 1
        zzzz(:,jmu) = zzzz(:,imu)
        slam(jmu) = slam(imu)
      end do
    end if
    do imu = 1, norbitals_new
      sqlami = slam(imu)**(-0.25d0)
      zzzz(:,imu) = zzzz(:,imu)*sqlami
    end do
    call dgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz, norbitals, zzzz, norbitals, a0, xxxx, norbitals)
    if (iqout .eq. 3) then
      do imu=1, norbitals
        xxxx(imu,:)=xxxx(imu,:)*ww(imu)
      end do
    endif
    do inu = 1, norbitals
      do imu = 1, norbitals
        sm12_save(imu,inu) = xxxx(imu,inu)
      end do
    end do
  else ! (if Kscf .eq. 1 .and iqout .ne. 3)
    xxxx(:,:) = sm12_save(:,:)
  end if
  if (iqout .ne. 3) then
    call dsymm ( 'R', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
    call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, zzzz, norbitals, a0, yyyy, norbitals )
  else
    call dgemm ( 'C', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
    call dgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, zzzz, norbitals, xxxx, norbitals, a0, yyyy, norbitals )
  endif
  lwork = 1
  deallocate (work)
  allocate (work(lwork))
  call dsyev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, -1, info)
  lwork = work(1)
  deallocate (work)
  allocate(work(lwork))
  call dsyev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, lwork, info )
  if (info .ne. 0) call diag_error (info, 0)
  eigen_k(1:norbitals,ikpoint) = eigen(:)
  if (iqout .ne. 2) blowre(:,:,ikpoint) = real(yyyy(:,:))
  if (iqout .ne. 2 .and. icluster .ne. 1) blowim(:,:,ikpoint) = 0
  call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  if (iqout .ne. 3) then
    call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  else
    call dgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  end if
  bbnkre(:,:,ikpoint) = real(zzzz(:,:))
  if (icluster .ne. 1) bbnkim(:,:,ikpoint) = 0
  deallocate (xxxx)
  deallocate (yyyy)
  deallocate (zzzz)
  if (iqout .eq. 3) then
    deallocate (ww)
    deallocate (ssss)
  endif
  deallocate (work)
  deallocate (iwork)
  write(*,*) 'XXX kspace_withno kbbnkre =' ,bbnkre
  return
end 

