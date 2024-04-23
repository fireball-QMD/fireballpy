subroutine louie ( x_try, x_old, beta, r2, iter, max_order, nmsh )
  use M_system
  implicit none
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real*8, intent(in) :: beta         ! Mixing factor - not used
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate? - not used
  real*8, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real*8, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real*8 new vector on output
  real*8, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2
  real*8, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged
  real*8 wgt0
  real*8 wgrad
  real*8 wgt1
  real*8, allocatable, dimension(:) :: deltaF
  real*8, allocatable, dimension(:) :: deltaX
  real*8, allocatable, dimension(:)   :: auxvec
  real*8, allocatable, dimension(:,:) :: amat
  real*8, allocatable, dimension(:,:) :: bmat
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
  write(*,*)'  '
  write(*,*)' Welcome to louie; mix charges to self-consistency.'
  wgt0 = 0.05d0
  wgrad = 0.01d0
  wgt1 = 1.0d0
  if(.not. allocated(Fv))then
     ! Why 2? We need the current and the one before.
     ! I guess we could do with just one, but it would make the code less legible.
     allocate (Fv(nmsh,2))
     allocate (Xv(nmsh,2))
     allocate (r2_sav(max_scf_iterations))
     allocate (betaInvH(nmsh,nmsh))
     allocate (gamaH(nmsh,nmsh))
  end if
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in Louie '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        stop
     end if
     Xv(:,1) = Xv(:,2)
     Fv(:,1) = Fv(:,2)
     Xv(:,2) = x_old(:)
     Fv(:,2) = x_try(:) - x_old(:)
     r2 = dot_product(Fv(:,2),Fv(:,2))
     r2 = r2 / nmsh
     mix_order = min(iter,max_order)
     if(r2 .lt. tr2) then
        x_old(:) = x_old(:) + beta*Fv(:,2)
        return
     end if
     if((mix_order .eq. 1)) then ! .or. (itipo .eq. 2 .and. mix_order .le. 3)
        x_old(:) = x_old(:) + beta*Fv(:,2)
        betaInvH(:,:) = 0.0d0
        gamaH(:,:) = 0.0d0
        aux = 1.0d0 / wgrad**2
        aux2 = wgrad**2 / wgt0
        do i = 1, nmsh
           betaInvH(i,i) = aux
           gamaH(i,i) = aux2
        end do
        return
     end if !mix_order .eq. 1
     allocate (deltaF(nmsh))
     allocate (deltaX(nmsh))
     allocate (auxvec(nmsh))
     allocate (amat(nmsh,nmsh))
     allocate (bmat(nmsh,nmsh))
     allocate (ipiv(nmsh))
     deltaX = Xv(:,2) - Xv(:,1)
     renorm = sqrt(dot_product( deltaX, deltaX ))
     deltaF = (Fv(:,2) - Fv(:,1)) / renorm !Eq A8
     deltaX = deltaX / renorm                      !Eq A7
        aux = 1.0d0
        aux2 = wgt1**2
        do i = 1,nmsh
           do j = 1, nmsh
              gamaH(i,j) = gamaH(i,j) -  aux2 * deltaF(i) * deltaX(j)  !Eq A15
              amat(i,j) = aux2 * deltaX(i) * deltaX(j)                 !Eq A16 first part
              aux = aux + aux2 * deltaX(i) * betaInvH(i,j) * deltaX(j) !Eq A16 second part
           end do
        end do
        bmat(:,:) = 0.0d0
        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,amat,nmsh,betaInvH,nmsh,0.0d0,bmat,nmsh)
        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,betaInvH,nmsh,bmat,nmsh,0.0d0,amat,nmsh)
        betaInvH(:,:) = betaInvH(:,:) - (1/aux) * amat(:,:)  !Eq A16 finish
        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,gamaH,nmsh,betaInvH,nmsh,0.0d0,amat,nmsh) !Eq A13, efficiently
     auxvec(:) = Fv(:,2)
     info = 0
     nrhs = 1
     lda = nmsh
     ldb = nmsh
     call dgesv(nmsh,nrhs,amat,lda,ipiv,auxvec,ldb,info) !get x**(m+1) - x**m
     x_old(:) = x_old(:) + auxvec(:)
     deallocate(deltaF,deltaX,auxvec)
     deallocate(amat,bmat)
     deallocate(ipiv)
     return
300  format (2x, ' norm =  ', f20.4)
305  format (2x, ' norm of difference between auxvecs =  ', f20.4)
301  format (2x, ' norm of the old vector =  ', f20.4)
302  format (2x, ' norm of the new guess =  ', f20.4)
308  format (2x, ' norm of difference between betaInvH, betaH^-1 =  ', f20.4)
242  format (2x, ' Selected method is:  ', a17)
   end subroutine louie

