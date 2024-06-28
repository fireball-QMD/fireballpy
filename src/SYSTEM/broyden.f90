subroutine broyden ( x_try, x_old, beta, r2, iter, max_order, nmsh)
  use M_system
  implicit none
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real*8, intent(in) :: beta         ! Mixing factor - used just for starting guess and simple mixing
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate? - not used
  real*8, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real*8, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real*8 new vector on output
  real*8, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2
  real*8, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged
  real*8, dimension(nmsh)   :: auxvec
  real*8, dimension(nmsh,nmsh) :: amat
  real*8 renorm
  real*8 aux
  real*8 aux2
  real*8 norm
  integer i,j,k
  integer mix_order ! Actual order used min(iter,max_order)
  integer, dimension(nmsh) :: ipiv
  integer info
  integer nrhs
  integer lda
  integer ldb
  integer ldx
  if(iter .eq. 1)then
     allocate (Fv(nmsh,2))
     allocate (Xv(nmsh,2))
     allocate (delX(nmsh,1), delF(nmsh,1))
     allocate (r2_sav(max_scf_iterations))
     allocate (RJac(nmsh,nmsh))
  end if
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in Broyden '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        deallocate(Fv, Xv, delF, delX, r2_sav, RJac)
        stop
     end if
     Xv(:,1) = Xv(:,2)
     Fv(:,1) = Fv(:,2)
     Xv(:,2) = x_old(:)
     Fv(:,2) = x_try(:) - x_old(:)
     r2 = dot_product(Fv(:,2), Fv(:,2))
     r2 = r2 / nmsh
     mix_order = min(iter,max_order)
     if((mix_order .eq. 1)) then
        x_old(:) = x_old(:) + beta*Fv(:,2)
        RJac(:,:) = 0.0d0
        aux = 1/beta
        do i = 1, nmsh
           RJac(i,i) = aux
        end do
        return
     end if ! mix_order .eq. 1
    if(r2 .lt. sigmatol) then
      scf_achieved = .true.
      deallocate(Fv, Xv, delF, delX, r2_sav, RJac)
      return
    end if
     delX(:,1) = Xv(:,2) - Xv(:,1)
     renorm = sqrt(dot_product( delX(:,1), delX(:,1) ))
     delF(:,1) = (Fv(:,2) - Fv(:,1)) / renorm !Eq A8
     delX(:,1) = delX(:,1) / renorm                      !Eq A7
     auxvec(:) = 0.0d0
     call dgemv('n',nmsh,nmsh,1.0d0,RJac,nmsh,delX(:,1),1,0.0d0,auxvec,1)
     do i = 1,nmsh
        do j = 1, nmsh
           RJac(i,j) = RJac(i,j) - ( delF(i,1) + auxvec(i) ) * delX(j,1) !Eq A10
        end do
     end do
     amat(:,:) = RJac(:,:)
     auxvec(:) = Fv(:,2)
     ipiv(:) = 0
     info = 0
     nrhs = 1
     lda = nmsh
     ldb = nmsh
     call dgesv(nmsh,nrhs,amat,lda,ipiv,auxvec,ldb,info) !get x^(m+1) - x^(m) from Eq A6
     x_old(:) = x_old(:) + auxvec(:)
   end subroutine broyden

