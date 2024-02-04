subroutine pulay ( x_try, x_old, beta, r2, iter, max_order, nmsh )
  use M_system
  implicit none
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real, intent(in) :: beta         ! Mixing factor
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate?
  real, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real new vector on output
  real, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2
  real, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged
  real, parameter :: A0 = 0.2d0  ! 0.2 can sometimes work better
  real, parameter :: q0 = 1.5d0  ! in angstrom^-1
  real, allocatable, dimension(:) :: deltaF
  real, allocatable, dimension(:) :: deltaX
  real, allocatable, dimension(:) :: metric   !for Eq.64
  real, allocatable, dimension(:) :: mixcoeff !alphas in Eq.52
  real, allocatable, dimension(:)     :: auxvec
  real, allocatable, dimension(:)     :: auxvec2
  real, allocatable, dimension(:,:) :: amat
  real, allocatable, dimension(:,:) :: bmat
  real, allocatable, dimension(:,:) :: cmat
  real, allocatable, dimension(:,:) :: idmat
  real renorm
  real aux
  real aux2
  real norm
  integer i,j,k
  integer mix_order ! Actual order used min(iter,max_order)
  integer, allocatable, dimension(:) :: ipiv
  integer info
  integer nrhs
  integer lda
  integer ldb
  integer ldx
  write(*,*)'  '
  write(*,*)' Welcome to pulay; mix charges to self-consistency.'
  if(.not. allocated(Fv))then
     allocate (Fv(nmsh,max_scf_iterations))
     allocate (Xv(nmsh,max_scf_iterations))
     allocate (r2_sav(max_scf_iterations))
  end if
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in pulay '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        stop
     end if
     Xv(:,iter) = x_old(:)
 !residual vector
     Fv(:,iter) = (x_try(:) - x_old(:))*mwe(:)        !residual vector
     r2 = dot_product(Fv(:,iter),Fv(:,iter))
     r2 = r2 / nmsh
     do i = 1,nmsh
       write (1000+iter,*) Fv(i,iter), drwe(i)
       write (2000+iter,*) Fv(i,iter), mwe(i)
     enddo
     mix_order = min(iter,max_order)
     if(r2 .lt. tr2) then
        x_old(:) = x_old(:) + beta*Fv(:,iter)/mwe(:)
        return
     end if
     if(mix_order .lt. 4) then
         write (*,*) ' Doing simple mixing', beta
         x_old(:) = (1.0d0-beta)*x_old(:) + beta*x_try(:)
        return
     end if
     allocate (mixcoeff(mix_order + 1)) !the last coefficient stands for the lagrangian multiplier for constraint of Eq.55
     allocate (amat(mix_order + 1, mix_order + 1))
     allocate(ipiv(mix_order + 1))
     amat = 0.0d0
     do i=1,mix_order
        do j=1,mix_order
           do k=1,nmsh
              amat(i,j) = amat(i,j) +  Fv(k,iter - mix_order + i) * Fv(k,iter - mix_order + j)
           end do
        end do
     end do
     amat(mix_order + 1,:) = -1.0d0
     amat(:,mix_order + 1) = -1.0d0
     amat(mix_order + 1, mix_order + 1) = 0.0d0
     mixcoeff(:) = 0.0d0
     mixcoeff(mix_order + 1) = -1.0d0
     ipiv = 0
     nrhs = 1
     lda = mix_order + 1
     ldb = mix_order + 1
     info = 0
     call dgesv(mix_order + 1,nrhs,amat,lda,ipiv,mixcoeff,ldb,info)
     write (*,*) 'sum_alpha =', sum(mixcoeff(1:mix_order))
     x_old(:) = 0.0d0
     do i = 1, mix_order
        x_old(:) = x_old(:) + mixcoeff(i) * Xv(:,i) !Eq.52
     end do
     deallocate(mixcoeff)
     deallocate(amat,ipiv)
     return
300  format (2x, ' norm =  ', f20.4)
305  format (2x, ' norm of difference between auxvecs =  ', f20.4)
301  format (2x, ' norm of the old vector =  ', f20.4)
302  format (2x, ' norm of the new guess =  ', f20.4)
308  format (2x, ' norm of difference between betaInvH, betaH^-1 =  ', f20.4)
242  format (2x, ' Selected method is:  ', a17)
   contains
     real function vecnorm(vec)
       real, intent(in), dimension(nmsh) :: vec
       integer i,j
       real r
       r = 0.0d0
       do i = 1,nmsh
          r = r + vec(i)**2
       end do
       r = sqrt(r)
       vecnorm = r
     end function  vecnorm
     real function maxnorm(mat)
       real, intent(in), dimension(nmsh,nmsh) :: mat
       integer i,j
       real m
       real r
       r = 0.0d0
       m = 0.0d0
       do i = 1,nmsh
          do j = 1,nmsh
             r = r + abs(mat(i,j))
          end do
          if (r .gt. m)  m = r
          r = 0.0d0
       end do
       maxnorm = m
     end function  maxnorm
end subroutine pulay

