subroutine pulay ( x_try, x_old, beta, r2, iter, max_order, nmsh )
  use M_system
  implicit none
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real*8, intent(in) :: beta         ! Mixing factor
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate?
  real*8, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real*8, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real*8 new vector on output
  real*8, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2
  real*8, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged
  real*8, parameter :: A0 = 0.2d0  ! 0.2 can sometimes work better
  real*8, parameter :: q0 = 1.5d0  ! in angstrom^-1
  real*8, allocatable, dimension(:) :: mixcoeff !alphas in Eq.52
  real*8, allocatable, dimension(:,:) :: amat
  real*8 renorm
  real*8 aux
  real*8 aux2
  real*8 norm
  integer i,j,k
  integer mix_order ! Actual order used min(iter,max_order)
  integer, allocatable, dimension(:) :: ipiv
  integer info
  integer nrhs
  integer lda
  integer ldb
  integer ldx
  if(iter .eq. 1) then
     allocate (Fv(nmsh,max_scf_iterations))
     allocate (Xv(nmsh,max_scf_iterations))
     allocate (r2_sav(max_scf_iterations))
  end if
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in pulay '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
       deallocate(Fv,Xv,r2_sav)
        stop
     end if
     Xv(:,iter) = x_old(:)
 !residual vector
     Fv(:,iter) = (x_try(:) - x_old(:))
     r2 = dot_product(Fv(:,iter),Fv(:,iter))
     r2 = r2 / nmsh
     mix_order = min(iter,max_order)
     print*, iter, abs(r2 - sigmatol)
     if(r2 .lt. sigmatol) then
       scf_achieved = .true.
       deallocate(Fv,Xv,r2_sav)
        return
     end if
     if(mix_order .lt. 4) then
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
     x_old(:) = 0.0d0
     do i = 1, mix_order
        x_old(:) = x_old(:) + mixcoeff(i) * Xv(:,i) !Eq.52
     end do
     deallocate(mixcoeff)
     deallocate(amat,ipiv)
     return
end subroutine pulay

