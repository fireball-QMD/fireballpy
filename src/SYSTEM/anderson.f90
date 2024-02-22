subroutine anderson ( x_try, x_old, beta, r2, iter, max_order, nmsh)
   use M_system
   implicit none
   integer, intent(in) :: nmsh      ! Size of vectors being optimized
   real, intent(in) :: beta         ! Mixing factor
   integer, intent(in) :: iter      ! iteration number
   integer, intent(in) :: max_order ! How far back do we go to extrapolate?
   real, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
   real, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real new vector on output
   real, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2 sigma
   real, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged
   real, allocatable, dimension(:,:) :: a_matrix  ! Eq. 5.17
   real, allocatable, dimension(:) :: delF_F      ! <delF|F> in Eq. 5.31
   real, allocatable, dimension(:) :: contribution
   integer iloop
   integer jloop
   integer mix_order ! Actual order used min(iter,max_order)
   integer lwork
   real, allocatable, dimension(:) :: work
   integer, allocatable, dimension(:) :: ipiv
   integer info
   if(.not. allocated(Fv))then
      allocate (Fv(nmsh,max_scf_iterations))
      allocate (Xv(nmsh,max_scf_iterations))
      allocate (delF(nmsh,max_scf_iterations))
      allocate (delX(nmsh,max_scf_iterations))
      allocate (r2_sav(max_scf_iterations))
   end if
   if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
      write (*,*) ' Stop in Anderson '
      write (*,*) ' max_order=', max_order
      write (*,*) ' iter=',iter
      write (*,*) ' max_scf_iterations=',max_scf_iterations
      stop
   end if
   Xv(:,iter) = x_old(:)
   Fv(:,iter) = x_try(:) - x_old(:)
   r2 = dot_product(Fv(:,iter),Fv(:,iter))
   r2 = r2 / nmsh
   mix_order = min(iter,max_order)
   if(mix_order .eq. 1 .or. r2 .lt. tr2) then
      x_old(:)=x_old(:) + beta*Fv(:,iter)
      return
   end if
   r2_sav(iter) = r2
   delF(:,iter-1) = Fv(:,iter) - Fv(:,iter-1) ! Eq. 5.6
   delX(:,iter-1) = Xv(:,iter) - Xv(:,iter-1) ! Eq. 5.5
   if (iter .gt. max_order .and. max_order .ge. 6) then
      if (r2_sav(iter-max_order) .lt. minval(r2_sav(iter-max_order+1:iter))) then
         r2_sav(iter-max_order+1) = r2_sav(iter-max_order)
         Fv(:,iter-max_order+1) = Fv(:,iter-max_order)
         Xv(:,iter-max_order+1) = Xv(:,iter-max_order)
         delX(:,iter-max_order+1) = Xv(:,iter-max_order+2) - Xv(:,iter-max_order+1)
         delF(:,iter-max_order+1) = Fv(:,iter-max_order+2) - Fv(:,iter-max_order+1)
      end if
   end if
888 allocate(a_matrix(iter-mix_order+1:iter-1,iter-mix_order+1:iter-1))
   do iloop = iter-mix_order+1, iter-1
      do jloop = iter-mix_order+1, iter-1
         a_matrix(iloop,jloop) = dot_product(delF(:,iloop),delF(:,jloop))
      end do
   end do
   allocate (delF_F(iter-mix_order+1:iter-1))
   do iloop = iter-mix_order+1, iter-1
      delF_F(iloop) = dot_product(delF(:,iloop),Fv(:,iter))  
   end do
   lwork = (mix_order-1)**2
   allocate (work(lwork))
   allocate (ipiv(mix_order-1))
   info=0
   call dsysv('U', mix_order-1, 1, a_matrix, mix_order-1, ipiv, delF_F, mix_order-1, work, lwork, info )
   if(info .ne. 0) then
      write (*,*) ' Error in Anderson, info =',info
      if(mix_order .le. 2)stop ! if you can't solve a 2x2 something is wrong
      mix_order=mix_order-1
      deallocate (work,ipiv,delF_F,a_matrix)
      goto 888 ! Try again with lower order
   end if
   x_old(:) = x_old(:) + beta*Fv(:,iter)  ! First-order term
   do iloop = iter-mix_order+1, iter-1
      x_old(:)=x_old(:) - delF_F(iloop)*(delX(:,iloop) + beta*delF(:,iloop))
   end do
   deallocate(delF_F,a_matrix,work,ipiv)
   return
 end subroutine anderson

