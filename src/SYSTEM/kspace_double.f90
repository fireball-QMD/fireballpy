subroutine kspace_double (ikpoint, sks)
  use M_system
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer, intent (in) :: ikpoint
  real(8), intent (in), dimension (3) :: sks
  real(8), parameter :: overtol = 1.0d-4
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
  real(8) dot
  real(8) sqlami
  real(8), dimension (norbitals) :: eigen
  real(8), dimension (norbitals) :: slam
  real(8), dimension (3) :: vec
  complex(16) a0
  complex(16) a1
  complex(16) phase
  complex(16), dimension (:, :), allocatable :: xxxx
  complex(16), dimension (:, :), allocatable :: yyyy
  complex(16), dimension (:, :), allocatable :: zzzz
  complex(16), dimension (:, :, :), allocatable, save :: sm12_save
  complex(16), dimension (:, :), allocatable :: ssss
  real(8), dimension (:), allocatable :: ww
  complex(16), allocatable, dimension (:) :: work
  real(8), allocatable, dimension (:) :: rwork
  integer, allocatable, dimension (:) :: iwork
  integer lwork, lrwork, liwork
  real(8) diff, imcoef, recoef
  a0 = cmplx(0.0d0,0.0d0)
  a1 = cmplx(1.0d0,0.0d0)
  ishort = 1
  allocate (xxxx(norbitals,norbitals))
  allocate (yyyy(norbitals,norbitals))
  allocate (zzzz(norbitals,norbitals))
  if (iqout .eq. 3) then
   allocate (ssss(norbitals,norbitals))
   allocate (ww(norbitals))
  endif

  lwork = 1
  allocate (work(lwork))
  lrwork = 3*norbitals - 2
  allocate (rwork(lrwork))

  if (.not. allocated(sm12_save))  allocate (sm12_save(norbitals,norbitals,nkpoints))
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
    vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
    dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
    phase = cmplx(cos(dot),sin(dot))
    do inu = 1, num_orb(in2)
     jnu = inu + degelec(jatom)
     do imu = 1, num_orb(in1)
      jmu = imu + degelec(iatom)
      zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
      yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
     end do !inu
    end do !imu
   end do !ineigh

   do ineigh = 1, neighPPn(iatom)
    mbeta = neighPP_b(ineigh,iatom)
    jatom = neighPP_j(ineigh,iatom)
    in2 = imass(jatom)
    vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
    dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
    phase = cmplx(cos(dot),sin(dot))
    do inu = 1, num_orb(in2)
     jnu = inu + degelec(jatom)
     do imu = 1, num_orb(in1)
      jmu = imu + degelec(iatom)
      yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*vnl(imu,inu,ineigh,iatom)
     end do !imu
    end do ! inu
   end do !inegh

  end do ! iatom                


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
  endif  

  if (Kscf .eq. 1 .or. iqout .eq. 3) then
   call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work, -1, rwork, info)
   lwork = work(1)
   deallocate (work)
   allocate (work(lwork))
   call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work, lwork, rwork , info)
   if (info .ne. 0) call diag_error (info, 0)
   mineig = 0
   do imu = 1, norbitals
    if (slam(imu) .lt. overtol) mineig = imu
   end do
   mineig = mineig + 1
   norbitals_new = norbitals + 1 - mineig
   if (norbitals_new .ne. norbitals) then
    write (*,*) ' Linear dependence encountered in basis set. '
    write (*,*) ' An overlap eigenvalue is very small. '
    write (*,*) norbitals - norbitals_new, ' vectors removed. '
    write (*,*) ' Spurious orbital energies near zero will '
    write (*,*) ' appear as a result of dropping these orbitals'
    write (*,*) ' You can change this by adjusting overtol in '
    write (*,*) ' kspace.f '
    write (*,*) '  '
    if(ishort .eq. 1) then    ! Don't print out again if done above
     write (*,*) '      The overlap eigenvalues: '
     write (*,*) ' ********************************************** '
     write (*,*) (slam(imu), imu = 1, norbitals)
    else          ! They asked for extra printout
     write(*,*) ' '
     write(*,*) ' Eigenvectors that correspond to eigenvalues'
     write(*,*) ' that were eliminated.  These might provide'
     write(*,*) ' insight into what atomic orbitals are causing'
     write(*,*) ' the problem.'
     write(*,*) ' '
     do imu = 1, mineig - 1
      write(*,*) ' eigenvector',imu
      do jmu = 1, norbitals
       write(*,*) jmu,' ',zzzz(jmu,imu)
      end do
     end do
    end if
    write (*,*) ' '
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
   call zgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz, norbitals, zzzz, norbitals, a0, xxxx, norbitals)
   if (iqout .eq. 3) then
    do imu=1, norbitals
     xxxx(imu,:)=xxxx(imu,:)*ww(imu)
    end do
   endif
   do inu = 1, norbitals
    do imu = 1, norbitals
     sm12_save(imu,inu,ikpoint) = xxxx(imu,inu)
    end do
   end do
  else ! (if Kscf .eq. 1 .and iqout .ne. 3)
   xxxx(:,:) = sm12_save(:,:,ikpoint)
  end if
  if (iqout .ne. 3) then
   call zhemm ( 'R', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
   call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, zzzz, norbitals, a0, yyyy, norbitals )
  else
   call zgemm ('C', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, ssss, norbitals, a0, zzzz, norbitals)
   call zgemm ('N', 'N', norbitals, norbitals, norbitals, a1, zzzz, norbitals, xxxx, norbitals, a0, zzzz, norbitals)
   call zgemm ( 'C', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
   call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, zzzz, norbitals, xxxx, norbitals, a0, yyyy, norbitals )
  endif     

  lwork = 1
  deallocate (work)
  allocate (work(lwork))
  call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, -1, rwork, info)
  lwork = work(1)
  deallocate (work)
  allocate (work(lwork))
  call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, lwork, rwork, info)
  if (info .ne. 0) call diag_error (info, 0)
  eigen_k(1:norbitals,ikpoint) = eigen(:)
  if (iqout .ne. 2) blowre(:,:,ikpoint) = real(yyyy(:,:), 8)
  if (iqout .ne. 2 .and. icluster .ne. 1) blowim(:,:,ikpoint) = aimag(yyyy(:,:))
  if (iqout .ne. 3) then
   call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  else
   call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  end if
  bbnkre(:,:,ikpoint) = real(zzzz(:,:), 8)
  if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(zzzz(:,:))
  deallocate (xxxx)
  deallocate (yyyy)
  deallocate (zzzz)
  if (iqout .eq. 3) then
    deallocate (ww)
    deallocate (ssss)
  endif
  deallocate (rwork)
  deallocate (work)
  return
end subroutine kspace_double

