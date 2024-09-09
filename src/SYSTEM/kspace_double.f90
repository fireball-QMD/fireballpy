subroutine kspace_double (ikpoint, sks)
  use iso_c_binding
  use M_system, only: iqout, icluster, igamma, natoms, ratom, degelec, imass, getlssh, getissh, Kscf, blowre, bbnkre, blowim, bbnkim, &
    & sm12_complex, eigen_k, norbitals, norbitals_new, getiatom, neigh_b, neigh_j, neighn, neighPPn, neighPP_b, neighPP_j, vnl, s_mat, &
    & h_mat, xl
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer(c_long), intent (in) :: ikpoint
  real(c_double), intent (in), dimension (3) :: sks
  real(c_double), parameter :: overtol = 1.0d-4
  integer(c_long) iatom
  integer(c_long) imu
  integer info
  integer(c_long) inu
  integer(c_long) in1, in2
  integer(c_long) ineigh
  integer(c_long) ishort
  integer(c_long) jatom
  integer(c_long) jmu
  integer(c_long) jnu
  integer(c_long) mbeta
  integer(c_long) mineig
  integer(c_long) lm
  real(c_double) dot
  real(c_double) sqlami
  real(c_double), dimension (norbitals) :: eigen
  real(c_double), dimension (norbitals) :: slam
  real(c_double), dimension (3) :: vec
  complex(c_double_complex) a0
  complex(c_double_complex) a1
  complex(c_double_complex) phase
  complex(c_double_complex), dimension (:, :), allocatable :: xxxx
  complex(c_double_complex), dimension (:, :), allocatable :: yyyy
  complex(c_double_complex), dimension (:, :), allocatable :: zzzz
  complex(c_double_complex), dimension (:, :), allocatable :: ssss
  real(c_double), dimension (:), allocatable :: ww
  complex(c_double_complex), allocatable, dimension (:) :: work
  real(c_double), allocatable, dimension (:) :: rwork
  integer(c_long) lwork, lrwork

  a0 = cmplx(0.0d0,0.0d0,c_double_complex)
  a1 = cmplx(1.0d0,0.0d0,c_double_complex)
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
    phase = cmplx(cos(dot),sin(dot),c_double_complex)
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
    phase = cmplx(cos(dot),sin(dot),c_double_complex)
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
   call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work, -1_c_long, rwork, info)
   lwork = nint(real(work(1), c_double), c_long)
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
     sm12_complex(imu,inu,ikpoint) = xxxx(imu,inu)
    end do
   end do
  else ! (if Kscf .eq. 1 .and iqout .ne. 3)
   xxxx(:,:) = sm12_complex(:,:,ikpoint)
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
  call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, -1_c_long, rwork, info)
  lwork = nint(real(work(1), c_double), c_long)
  deallocate (work)
  allocate (work(lwork))
  call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, lwork, rwork, info)
  if (info .ne. 0) call diag_error (info, 0)
  eigen_k(1:norbitals,ikpoint) = eigen(:)
  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    blowre(:,:,ikpoint) = real(yyyy(:,:), c_double)
    if (igamma .eq. 0) then
      blowim(:,:,ikpoint) = aimag(yyyy(:,:))
    end if
  end if
  if (iqout .ne. 3) then
   call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  else
   call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx, norbitals, yyyy, norbitals, a0, zzzz, norbitals )
  end if
  bbnkre(:,:,ikpoint) = real(zzzz(:,:), c_double)
  if (icluster .eq. 0 .and. igamma .eq. 0) bbnkim(:,:,ikpoint) = aimag(zzzz(:,:))
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
