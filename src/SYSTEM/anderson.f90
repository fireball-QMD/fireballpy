subroutine anderson ( x_try, x_old, nmsh )
  use iso_c_binding
  use M_system, only: ialgmix, idmix, max_scf_iterations, bmix, sigma, sigmatol, sigmabest, w02, scf_achieved, Kscf, Fv, Xv, delF, delX, &
    & r2_sav, wi, x_best
  implicit none
  integer(c_long), intent(in) :: nmsh    ! Size of vectors being optimized
  real(c_double), intent(in), dimension(nmsh) :: x_try ! potential new vector on input
  real(c_double), intent(inout), dimension(nmsh) :: x_old ! old vector in input, real(c_double) new vector on output
  integer(c_long) :: iloop, jloop, mix_order, imix, Kscfm1
  real(c_double), allocatable, dimension(:) :: delF_F    ! <delF|F> in Eq. 5.31
  real(c_double), allocatable, dimension(:,:) :: a_matrix
  ! BLAS stuff
  integer(c_long) :: lwork
  integer :: info
  integer(c_long), allocatable, dimension(:) :: ipiv
  real(c_double), allocatable, dimension(:) :: work

  if(Kscf .eq. 1)then
    ! wi need to be allocated somewhere else
    allocate (Fv(nmsh,max_scf_iterations+1))
    allocate (Xv(nmsh,max_scf_iterations+1))
    allocate (delF(nmsh,max_scf_iterations))
    allocate (delX(nmsh,max_scf_iterations))
    allocate (r2_sav(max_scf_iterations))
    allocate (x_best(nmsh))

    Xv(:,1) = x_old(:)
    Fv(:,1) = x_try(:) - x_old(:)
    sigma = dot_product(Fv(:,1), Fv(:,1)) / nmsh
    if(sigma .lt. sigmatol) then
      scf_achieved = .true.
      deallocate(Fv, Xv, delF, delX, r2_sav, x_best, wi)
      return
    end if
    x_old(:) = x_old(:) + bmix*Fv(:,1)
    r2_sav(1) = sigma
    sigmabest = sigma
    return
  end if

  Xv(:,Kscf) = x_old(:) ! last element is x_old, rest are |U>
  Fv(:,Kscf) = x_try(:) - x_old(:)
  sigma = dot_product(Fv(:,Kscf), Fv(:,Kscf)) / nmsh
  if(sigma .lt. sigmatol) then
    scf_achieved = .true.
    deallocate(Fv, Xv, delF, delX, r2_sav, wi, x_best)
    return
  else if(Kscf .eq. max_scf_iterations+1 .and. sigma .gt. sigmabest) then
    x_old(:) = x_best(:)
    deallocate(Fv, Xv, delF, delX, r2_sav, wi, x_best)
    return
  else if(sigma .lt. sigmabest) then
    x_best(:) = x_old(:)
    sigmabest = sigma
  end if
  r2_sav(Kscf) = sigma

  Kscfm1 = Kscf - 1
  if (ialgmix .eq. 2) then
    wi(Kscfm1) = 1.0d0 / r2_sav(Kscfm1)
  end if

  delX(:,Kscfm1) = wi(Kscfm1)*(x_old(:) - Xv(:,Kscfm1))
  delF(:,Kscfm1) = wi(Kscfm1)*(Fv(:,Kscf) - Fv(:,Kscfm1))
  Xv(:,Kscfm1) = bmix*delF(:,Kscfm1) + delX(:,Kscfm1)

  imix = max(1, Kscf - idmix + 1)
  mix_order = min(Kscf, idmix) - 1
  if (imix .gt. 1) then
    if (r2_sav(imix-1) .lt. minval(r2_sav(imix:Kscf))) then
      r2_sav(imix) = r2_sav(imix-1)
      Fv(:,imix) = Fv(:,imix-1)
      Xv(:,imix) = Xv(:,imix-1)
      delX(:,imix) = delX(:,imix-1)
      delF(:,imix) = delF(:,imix-1)
    end if
  end if

  allocate (delF_F(imix:Kscfm1))
  allocate (a_matrix(imix:Kscfm1,imix:Kscfm1))
  do jloop = imix,Kscfm1
    delF_F(jloop) = dot_product(delF(:,jloop), Fv(:,Kscf))
    do iloop = imix,jloop
      a_matrix(iloop,jloop) = dot_product(delF(:,iloop), delF(:,jloop))
    end do
    a_matrix(jloop,jloop) = a_matrix(jloop,jloop) + w02
  end do
  allocate(work(1), ipiv(mix_order))
  call dsysv('U',mix_order,1,a_matrix,mix_order,ipiv,delF_F,mix_order,work,-1_c_long,info)
  lwork = nint(work(1), c_long)
  deallocate(work)
  allocate(work(lwork))
  call dsysv('U',mix_order,1,a_matrix,mix_order,ipiv,delF_F,mix_order,work,lwork,info)

  x_old(:) = x_old(:) + bmix*Fv(:,Kscf)
  do jloop=imix,Kscfm1
    x_old = x_old(:) - delF_F(jloop)*Xv(:,jloop)
  end do
  deallocate (delF_F,a_matrix,work,ipiv)
end subroutine anderson
