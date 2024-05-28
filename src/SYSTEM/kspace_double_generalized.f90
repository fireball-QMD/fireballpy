subroutine kspace_double_generalized (ikpoint, sks)
  use M_system
  use M_fdata, only: num_orb, Qneutral
  implicit none
  integer, intent (in) :: ikpoint
  real*8, intent (in), dimension (3) :: sks
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
  real*8 dot
  real*8 sqlami
  real*8, dimension (norbitals) :: eigen
  real*8, dimension (norbitals) :: slam
  real*8, dimension (3) :: vec
  complex*16 a0
  complex*16 a1
  complex*16 phase
  complex*16, dimension (:, :), allocatable :: xxxx
  complex*16, dimension (:, :), allocatable :: yyyy
  complex*16, dimension (:, :), allocatable :: zzzz
  complex*16, dimension (:, :), allocatable :: ssss
  real*8, dimension (:), allocatable :: ww
  complex*16, allocatable, dimension (:) :: work
  real*8, allocatable, dimension (:) :: rwork
  integer, allocatable, dimension (:) :: iwork
  integer lwork, lrwork, liwork

  real*8 diff, imcoef, recoef
  a0 = cmplx(0.0d0,0.0d0)
  a1 = cmplx(1.0d0,0.0d0)
  ishort = 1
  allocate (xxxx(norbitals,norbitals))
  allocate (yyyy(norbitals,norbitals))
  allocate (zzzz(norbitals,norbitals))

  lwork = norbitals*norbitals ! Use xxxx, yyyy and zzzz for work area
  lrwork = 3*norbitals
  allocate (rwork(lrwork))
  allocate (work(lwork))

  zzzz = a0
  yyyy = a0
  xxxx = a0

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

  ! optionally, save Hamiltonian & Overlap for DOS
  xxxx = yyyy


  norbitals_new = norbitals


  ! xxxx = unused (used as complex workspace in cheev call below)
  ! zzzz = Overlap in AO basis
  ! yyyy = Hamiltonian in AO basis
  
  ! DIAGONALIZE THE HAMILTONIAN.
  ! ****************************************************************************
  !
  ! ZHEGV  -  compute all the eigenvalues, and optionally, the eigenvectors
  !   A*x=(lambda)*B*x
  ! Eigenvectors are needed to calculate the charges and for forces!
  ! Compute the eigenvalues and the eigenvectors of a  complex  generalized
  ! Hermitian-definite  eigenproblem,  of the form A*x=(lambda)*B*x,

  call zhegv (1, 'V', 'U', norbitals, yyyy, norbitals, zzzz,    &
           &  norbitals, eigen, work, lwork, rwork , info)

  if (info .ne. 0) call diag_error (info, 0)

  eigen_k(1:norbitals,ikpoint) = eigen(:)
  bbnkre(:,:,ikpoint) = real(yyyy(:,:))

  if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(yyyy(:,:))

  deallocate (xxxx)
  deallocate (yyyy)
  deallocate (zzzz)

  deallocate (rwork)
  deallocate (work)

  return
end subroutine 

