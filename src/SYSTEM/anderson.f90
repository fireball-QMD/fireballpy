subroutine anderson ( x_try, x_old, beta, r2, iter, max_order, nmsh)
  use M_system
  implicit none
  integer, intent(in) :: nmsh    ! Size of vectors being optimized
  real*8, intent(in) :: beta      ! Mixing factor
  integer, intent(in) :: iter    ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate?
  real*8, intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real*8, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real*8 new vector on output
  real*8, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2
  integer :: iloop, jloop, retry, mix_order, imix, iterm1, mixm1
  real*8, allocatable, dimension(:) :: delF_F    ! <delF|F> in Eq. 5.31
  real*8, allocatable, dimension(:,:) :: a_matrix

  !real*8, parameter :: w02 = 0.0001
  real*8, parameter :: w02 = 0.0
  ! BLAS stuff
  integer :: lwork, info
  integer, allocatable, dimension(:) :: ipiv
  real*8, allocatable, dimension(:) :: work

  if(iter .eq. 1)then
    allocate (Fv(nmsh,max_scf_iterations))
    allocate (Xv(nmsh,max_scf_iterations))
    allocate (delF(nmsh,max_scf_iterations))
    allocate (delX(nmsh,max_scf_iterations))
    allocate (r2_sav(max_scf_iterations))

    Xv(:,1) = x_old(:)
    Fv(:,1) = x_try(:) - x_old(:)
    r2 = dot_product(Fv(:,1), Fv(:,1)) / nmsh
    if(r2 .lt. sigmatol) then
      scf_achieved = .true.
      deallocate(Fv, Xv, delF, delX, r2_sav)
      return
    end if
    x_old(:) = x_old(:) + beta*Fv(:,1)
    r2_sav(1) = r2
    return
  end if

  ! Need to be changed
  if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
    write (*,*) ' Stop in Anderson '
    write (*,*) ' max_order=', max_order
    write (*,*) ' iter=',iter
    write (*,*) ' max_scf_iterations=',max_scf_iterations
    deallocate(Fv, Xv, delF, delX, r2_sav)
    stop
  end if

  Xv(:,iter) = x_old(:) ! last element is x_old, rest are |U>
  Fv(:,iter) = x_try(:) - x_old(:)
  r2 = dot_product(Fv(:,iter), Fv(:,iter)) / nmsh
  if(r2 .lt. sigmatol) then
    scf_achieved = .true.
    deallocate(Fv, Xv, delF, delX, r2_sav)
    return
  end if
  r2_sav(iter) = r2

  iterm1 = iter - 1
  delX(:,iterm1) = (x_old(:) - Xv(:,iterm1))!/(nmsh*r2_sav(iterm1))
  delF(:,iterm1) = (Fv(:,iter) - Fv(:,iterm1))!/(nmsh*r2_sav(iterm1))
  Xv(:,iterm1) = beta*delF(:,iterm1) + delX(:,iterm1)

  imix = max(1, iter - max_order + 1)
  mix_order = min(iter, max_order) - 1
  if (iter .gt. max_order .and. max_order .ge. 6) then
    if (r2_sav(imix-1) .lt. minval(r2_sav(imix:iter))) then
      r2_sav(imix) = r2_sav(imix-1)
      Fv(:,imix) = Fv(:,imix-1)
      Xv(:,imix) = Xv(:,imix-1)
      delX(:,imix) = delX(:,imix-1)
      delF(:,imix) = delF(:,imix-1)
    end if
  end if

  allocate (delF_F(imix:iterm1))
  allocate (a_matrix(imix:iterm1,imix:iterm1))
  do jloop = imix,iterm1
    delF_F(jloop) = dot_product(delF(:,jloop), Fv(:,iter))
    do iloop = imix,jloop
      a_matrix(iloop,jloop) = dot_product(delF(:,iloop), delF(:,jloop))
    end do
    a_matrix(jloop,jloop) = a_matrix(jloop,jloop) + w02
  end do
  allocate(work(1), ipiv(mix_order))
  call dsysv('U',mix_order,1,a_matrix,mix_order,ipiv,delF_F,mix_order,work,-1,info)
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))
  call dsysv('U',mix_order,1,a_matrix,mix_order,ipiv,delF_F,mix_order,work,lwork,info)

  x_old(:) = x_old(:) + beta*Fv(:,iter)
  do jloop=imix,iterm1
    x_old = x_old(:) - delF_F(jloop)*Xv(:,jloop)
  end do
  deallocate (delF_F,a_matrix,work,ipiv)
end subroutine anderson
