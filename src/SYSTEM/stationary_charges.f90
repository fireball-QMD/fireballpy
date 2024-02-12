subroutine STATIONARY_CHARGES()                   
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom                        
  integer ikpoint                    
  integer imu,inu                   
  integer in1, in2
  integer issh, jssh                         
  integer ineigh, jatom
  integer mbeta
  integer noccupy                    
  integer mqn
  real,dimension(nssh_tot,nssh_tot) :: A
  real,dimension(nssh_tot) :: c, SQ ! carga
  real,dimension(nssh_tot) :: LB, UB, nalpha
  real :: diff_err,Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real Ntot
  real auxgS
  integer :: beta, alpha, ialp, ina, matom
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
  do ialp = 1, natoms
    ina = imass(ialp)
    !gvhxc
    do issh = 1, nssh(ina)
      alpha = alpha + 1 ! transform to one index
      beta = 0
      do iatom = 1, natoms
        in1 = imass(iatom)
        matom = neigh_self(iatom)            
        inumorb = 1 ! counter for number of orbitals in atom iatom
        do issh1 = 1, nssh(in1)
          beta = beta + 1 ! transform to one index
          ! Spherical approximation to matrix elements:
          l = lssh(issh1,in1)
          auxgS = 0.0d0
          mu_min = inumorb
          mu_max = mu_min+2*l
          do imu = mu_min, mu_max 
            auxgS =  auxgS  +  gvhxc(imu,imu,issh,ialp,matom,iatom)
          end do ! end do imu = mu_min, mu_max 
          auxgS = auxgS/(2*l+1)  ! 4*pi??
          !M(alpha,beta) =  auxgS !gvhxcs(issh1,issh,iatom,ialp)          
          A(alpha,beta) = auxgS 
          inumorb = inumorb + 2*l+1
        end do ! end do issh1
        do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom) 
          do imu = 1, num_orb(in1)
            do inu = 1, num_orb(in2)
              c(alpha) = c(alpha) + rho(imu,inu,ineigh,iatom)*gvhxc(imu,inu,issh,ialp,ineigh,iatom)
            end do ! end do inu
          end do ! end do imu
        end do ! end do ineigh
      end do ! end do iatom
    end do ! end do issh
  end do ! end do ialp
  !SOLVE SYSTEM Mx = B.  x are the charges
  !M(nssh_tot+1,nssh_tot+1) = 0
  !B(1,nssh_tot+1) = ztot
  !do alpha = 1,nssh_tot
  ! M(nssh_tot+1,alpha) = 1
  ! M(alpha,nssh_tot+1) = 1
  !end do
  !LWMAX = 100
  !call ssysv( 'U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, lwork, info )
  !call sgesv(nssh_tot,1,M,nssh_tot,ipiv,B,nssh_tot,info )
  alpha = 0
  Ntot = 0.0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      alpha = alpha + 1
      SQ(alpha) = Qin(issh,iatom)
      LB(alpha) = 0.00
      if ( lssh(issh,in1) .eq. 0 ) then 
        UB(alpha) = 2.00
        nalpha(alpha) = 1.00 
      end if
      if ( lssh(issh,in1) .eq. 1 ) then
        UB(alpha) = 6.00
        nalpha(alpha) = 3.00 
      end if
      if ( lssh(issh,in1) .eq. 2 ) then
        UB(alpha) = 10.00
        nalpha(alpha) = 5.00 
      end if 
      !descomentar para anular UB y nalpha
      !UB(alpha) = 100.00
      !nalpha(alpha) = 1.00
      Ntot = Ntot + Qin(issh,iatom)       
    end do ! end do issh
  end do ! end do iatom
  call step_size(nssh_tot,A,c,SQ, LB, UB, nalpha) !,,LB,UB) B(1,alpha
  Ntot=0
  do alpha = 1, nssh_tot
    Ntot = Ntot + SQ(alpha)
  end do
  diff_err=Ep2(SQ,A,c,nssh_tot,nalpha)
  alpha = 0 
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      alpha = alpha + 1
      Qout(issh,iatom) = SQ(alpha) 
    end do ! end do issh
  end do ! end do iatom
end subroutine STATIONARY_CHARGES

subroutine step_size(nssh_tot,A,c,Q,LB,UB,nalpha) !,LB,UB,nstep)
  integer, intent(in) :: nssh_tot
  real,intent(in),dimension(nssh_tot,nssh_tot) :: A
  real,intent(in),dimension(nssh_tot) :: c
  real,intent(inout),dimension(nssh_tot) :: Q
  real,intent(in),dimension(nssh_tot) :: LB, UB, nalpha
  integer :: nstep = 0
  logical, dimension(nssh_tot) :: cero
  integer :: i, j, k,nceros
  real :: err0,err1,dp,g0,dqmax,x,y1,y2,y3,gmod,qmod
  real,dimension(nssh_tot) ::Q0,g,fc
  real,dimension(nssh_tot,nssh_tot) :: F
  real :: dq , diff_err
  logical :: g_err, igualceros
  !========== control ==========
  real :: tol_err=0.001, nstepmax=2000
  real :: tol = 0.001
  !=============================
  Q0=Q
  F=0
  do i=1,nssh_tot
      do j=1, nssh_tot
         do k=1, nssh_tot
              F(i,j)=F(i,j)+2*A(k,i)*nalpha(k)*A(k,j)
         end do
      end do
  end do
  fc=0
  do i=1,nssh_tot
      do j=1, nssh_tot
         fc(i)=fc(i)-2*c(j)*nalpha(j)*A(j,i)
      end do
  end do
  nstep=0
  diff_err=100.0
  do while (diff_err > tol_err ) ! .and. nstep .lt. nstepmax) 
      nstep=nstep+1
      ! proyectamos en plano Qtot = cte
      g=0.00
      do i=1,nssh_tot
         do j=1, nssh_tot
            g(i)=g(i)-F(i,j)*Q(j)
         end do
         g(i)=g(i)-fc(i)
      end do
      g0=0.00
      do i=1,nssh_tot
            g0=g0+g(i)
      end do
      g0=g0/nssh_tot
      do i=1,nssh_tot
         g(i)=g(i)-g0
      end do
      !proyectamos en las componentes q<0
      do i=1,nssh_tot
         if((Q0(i) < LB(i)+tol .and. g(i) < 0.00 ) .or. (Q0(i) > UB(i)-tol .and. g(i) > 0.00 )) then
            nceros=nceros+1
            g(i)=0.00
            cero(i)=.True.
         else
            cero(i)=.False.
            g0=g0+g(i)
         end if
      end do
      igualceros=.False.
      do while (igualceros .eqv. .False.)
         call getceros(nssh_tot,Q0,g,LB,UB,tol,cero,igualceros)
      end do
      nceros=0
      g0=0.00
      do i=1,nssh_tot
         if(cero(i) .eqv. .True.) then
            nceros=nceros+1
            g(i)=0.00
         else
            g0=g0+g(i)
         end if
      end do
      g0=g0/(nssh_tot-nceros)
      do i=1,nssh_tot
         if (cero(i) .eqv. .False.) then
            g(i)=g(i)-g0
         else
            g(i)=0.00
         end if
      end do
      gmod=0.00
      do i=1,nssh_tot
         gmod=gmod+g(i)**2
      end do
      gmod=gmod**0.5
      qmod=0.00
      do i=1,nssh_tot
         qmod=qmod+Q(i)**2
      end do
      qmod=qmod**0.5
      dq=qmod/gmod
      !obtenemos el min
      x=get_min_parabola(-dq,0,dq,Ep2(Q0-dq*g,A,c,nssh_tot,nalpha),Ep2(Q0,A,c,nssh_tot,nalpha),Ep2(Q0+dq*g,A,c,nssh_tot,nalpha))
      Q=Q0+x*g
      diff_err=Ep2(Q0,A,c,nssh_tot,nalpha)-Ep2(Q,A,c,nssh_tot,nalpha)
      dqmax=1.0E+10
      k=0
      do i=1,nssh_tot
         if(cero(i) .eqv. .False.) then  !ojo !!!!
            if(Q(i) < LB(i)) then
      if (dqmax > -Q0(i)/g(i)) then
         dqmax=-Q0(i)/g(i)
         k=k+1
      end if
            end if
            if(Q(i) > UB(i)) then
      if (dqmax > (UB(i)-Q0(i))/g(i)) then
         dqmax=(UB(i)-Q0(i))/g(i)
         k=k+1
      end if
            end if
         end if
      end do
      if(k .eq.0 ) then
         Q=Q
      else
         Q=Q0+dqmax*g
      end if
      Q0=Q
  end do
  aux=0
  do i = 1, nssh_tot
      aux=aux+Q(i)
  end do
end

real function Ep2(q,A,c,nssh_tot,nalpha)
  integer, intent(in) :: nssh_tot
  real,intent(in),dimension(nssh_tot,nssh_tot) :: A
  real,intent(in),dimension(nssh_tot) :: q,c,nalpha
  integer i,j,k
  real aux
  !Ep=(Aq-c)**2
  Ep2=0.00
  do i = 1, nssh_tot
      aux=0.00
      do j = 1, nssh_tot
         aux=aux+A(i,j)*q(j)
      end do
      aux=aux-c(i)
      aux=aux**2*nalpha(i)
      Ep2=Ep2+aux
  end do
end function Ep2


subroutine getceros(nssh_tot,Q0,g,LB,UB,tol,cero,igualceros)
  integer, intent(in) :: nssh_tot
  real, intent(in) :: tol
  real,intent(in),dimension(nssh_tot) :: Q0
  real,intent(in),dimension(nssh_tot) :: LB, UB
  real,intent(in),dimension(nssh_tot) :: g
  logical,intent(inout) ,dimension(nssh_tot) :: cero
  logical, intent(inout) :: igualceros
  real,dimension(nssh_tot) :: gaux
  real :: g0
  integer :: i, j,nceros
  real :: aux
  nceros=0
  g0=0.00
  gaux=g
  !write(*,'(A6,<nssh_tot>F12.3,A3)')'g0=(/ ',(gaux(i),i = 1, nssh_tot),' /)'
  !write(*,'(A6,<nssh_tot>F12.3,A3)')'Q0=(/ ',(Q0(i),i = 1, nssh_tot),' /)'
  do i=1,nssh_tot
      if(cero(i) .eqv. .True.) then
         nceros=nceros+1
         gaux(i)=0.00
      else
         g0=g0+g(i)
      end if
  end do
  g0=g0/(nssh_tot-nceros)
  do i=1,nssh_tot
      if (cero(i) .eqv. .False.) then
         gaux(i)=gaux(i)-g0
      else
         gaux(i)=0.00
      end if
  end do
  igualceros=.True.
  do i=1,nssh_tot
      if(cero(i)  .eqv. .False.)then
         if((Q0(i) < LB(i)+tol .and. gaux(i) < 0.00 ) .or. (Q0(i) > UB(i)-tol .and. gaux(i) > 0.00 )) then
            nceros=nceros+1
            cero(i)=.True.
            igualceros=.False.
         end if
      end if
  end do
end subroutine getceros
