subroutine denmat ()
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom
  integer iband
  integer ikpoint
  integer imu, inu
  integer ineigh
  integer in1, in2
  integer iorbital
  integer issh, jssh
  integer jatom
  integer jneigh
  integer mqn
  integer mbeta
  integer mmu
  integer noccupy
  integer nnu
  integer :: info, lwork
  integer, dimension(100) :: work
  real aux1, aux2, aux3
  real deltae
  real dot
  real gutr
  real pcharge
  real ztest
  real checksum
  real Wmu
  real, dimension (natoms) :: pqmu
  real, dimension (natoms) :: QoutTot
  real, dimension (3) :: vec
  complex ai
  complex phase, phasex
  complex step1, step2
  logical read_occupy
  ai = cmplx(0.0d0,1.0d0)
  rhoPP = 0.0d0
  !AQUI  inquire (file = 'OCCUPATION', exist = read_occupy)

  !Get the Fermi energy.
  call fermie ()
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      do ikpoint = 1, nkpoints
        dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) + special_k(3,ikpoint)*vec(3)
        phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
        if (icluster .ne. 1) then
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*(bbnkre(mmu,iband,ikpoint) - ai*bbnkim(mmu,iband,ikpoint))
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*(bbnkre(nnu,iband,ikpoint) + ai*bbnkim(nnu,iband,ikpoint))
                  gutr = real(step2)
                  rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                  cape(imu,inu,ineigh,iatom) = cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
                end do
              end do
            end if
          end do
        else
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*bbnkre(mmu,iband,ikpoint)
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*bbnkre(nnu,iband,ikpoint)
                  gutr = real(step2)
                  rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                  cape(imu,inu,ineigh,iatom) = cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do

  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighPPn(iatom)
      mbeta = neighPP_b(ineigh,iatom)
      jatom = neighPP_j(ineigh,iatom)
      in2 = imass(jatom)
      vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
      do ikpoint = 1, nkpoints
        dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)  + special_k(3,ikpoint)*vec(3)
        phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
        if (icluster .ne. 1) then
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*(bbnkre(mmu,iband,ikpoint) - ai*bbnkim(mmu,iband,ikpoint))
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*(bbnkre(nnu,iband,ikpoint)  + ai*bbnkim(nnu,iband,ikpoint))
                  gutr = real(step2)
                  rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom) + gutr
                end do
              end do
            end if
          end do
        else
          do iband = 1, norbitals_new
            if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*bbnkre(mmu,iband,ikpoint)
                do inu = 1, num_orb(in2)
                  nnu = inu + degelec(jatom)
                  step2 = step1*bbnkre(nnu,iband,ikpoint)
                  gutr = real(step2)
                  rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom)  + gutr
                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do

  if (iqout .eq. 1 .or. iqout .eq. 3) then
    Qout = 0.0d0
    QLowdin_TOT = 0.0d0
    call LOWDIN_CHARGES()
  end if !iqout = 1,3
  if (iqout .eq. 2) then
    Qout = 0.0d0
    QMulliken_TOT = 0.0d0
    call MULLIKEN_CHARGES()
  end if !iqout = 2
  if (iqout .eq. 4) then
    Qout = 0.0d0                        
    QMulliken_TOT = 0.0d0 
    call MULLIKEN_DIPOLE_CHARGES()
  end if !iqout = 4
  if (iqout .eq. 6) then
    call STATIONARY_CHARGES()
  end if !iqout = 6
  if (iqout .eq. 7) then
    Qout = 0.0d0                        
    QMulliken_TOT = 0.0d0 
    call MULLIKEN_DIPOLE_CHARGES()
    call Dipole_proyection()
    do iatom = 1, natoms
      in1 = imass(iatom)
      QoutTot(iatom) = 0.0d0
      do imu = 1,nssh(in1)
        QoutTot(iatom) = QoutTot(iatom)+Qout(imu,iatom)
      end do ! end do imu
    end do !end do iatom = 1,natoms
    do iatom = 1, natoms
      in1 = imass(iatom)
      do imu = 1,nssh(in1)
        Qout(imu,iatom) = (dq_DP(iatom)/QoutTot(iatom))*Qout(imu,iatom) + Qout(imu,iatom)
      end do ! end do imu
    end do !end do iatom = 1,natoms
  end if !iqout = 7
  ebs = 0.0d0
  ztest = 0.0d0
  do ikpoint = 1, nkpoints
    do iorbital = 1, norbitals_new
      if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
       ebs = ebs + weight_k(ikpoint)*spin*eigen_k(iorbital,ikpoint)*foccupy(iorbital,ikpoint)
       ztest = ztest + weight_k(ikpoint)*spin*foccupy(iorbital,ikpoint)
      end if
    end do
  end do
  if (abs(ztest - ztot) .gt. 1.0d-02) then
    write (*,*) ' *************** error *************** '
    write (*,*) ' ztest = ', ztest, ' ztot = ', ztot
    write (*,*) ' In denmat.f - ztest .ne. ztot! '
    stop
  end if
  return
end subroutine denmat

subroutine LOWDIN_CHARGES()
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom
  integer ikpoint
  integer imu, inu
  integer in1, in2
  integer issh, jssh,mmu
  integer noccupy
  integer mqn
  integer iorbital
  real aux1, aux2, aux3
  Qout = 0.0d0
  QLowdin_TOT = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ikpoint = 1, nkpoints
       aux1 = weight_k(ikpoint)*spin
       do iorbital = 1, norbitals
         if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
          aux2 = aux1*foccupy(iorbital,ikpoint)
          imu = 0
          do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
              imu = imu + 1
              mmu = imu + degelec(iatom)
              if (icluster .ne. 1) then
                aux3 = aux2*(blowre(mmu,iorbital,ikpoint)**2  + blowim(mmu,iorbital,ikpoint)**2)
              else
                aux3 = aux2*blowre(mmu,iorbital,ikpoint)**2
              end if
              Qout(issh,iatom) = Qout(issh,iatom) + aux3
              QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3
            end do
          end do
        end if
      end do !kpoints
    end do
  end do !atoms
end subroutine LOWDIN_CHARGES 

subroutine MULLIKEN_CHARGES()                
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom                        
  integer ikpoint                    
  integer imu, inu                   
  integer in1, in2                   
  integer issh, jssh
  integer ineigh ,jatom,jneigh                   
  integer noccupy    
  integer mqn                          
  real, dimension (numorb_max, natoms) :: QMulliken
  QMulliken = 0.0d0                         
  do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
         jatom = neigh_j(ineigh,iatom)
         in2 = imass(jatom)
         jneigh = neigh_back(iatom,ineigh)
         do imu = 1, num_orb(in1)
            do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)+ 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
            end do
         end do
      end do !ineig
      imu = 0
      do issh = 1, nssh(in1)
         do mqn = 1, 2*lssh(issh,in1) + 1
            imu = imu + 1
            Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
         end do
         QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(issh,iatom)
      end do
  end do    !atoms
end subroutine MULLIKEN_CHARGES 

subroutine MULLIKEN_DIPOLE_CHARGES()                   
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom                        
  integer ikpoint                    
  integer imu, inu                   
  integer in1, in2                   
  integer issh, jssh
  integer ineigh , jatom,jneigh                        
  integer noccupy    
  integer mqn                          
  real y
  real, dimension (numorb_max, natoms) :: QMulliken
  real, dimension (3) :: vec, r1, r2, r21
  QMulliken = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do ineigh = 1, neighn(iatom)
       jatom = neigh_j(ineigh,iatom)
       in2 = imass(jatom)
       r2(:) = ratom(:,jatom)
       ! Find r21 = vector pointing from r1 to r2, the two ends of the
       ! bondcharge, and the bc distance, y
       r21(:) = r2(:) - r1(:)
       y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
       jneigh = neigh_back(iatom,ineigh)
       do imu = 1, num_orb(in1)
          do inu = 1, num_orb(in2)
            QMulliken(imu,iatom) = QMulliken(imu,iatom) + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
          end do
       end do
       ! dipole correction. Only if the two atoms are different
       if (y .gt. 1.0d-05) then
          do imu = 1, num_orb(in1)
            do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)+ (-rho(imu,inu,ineigh,iatom)*dip(imu,inu,ineigh,iatom)+ rho(inu,imu,jneigh,jatom)*dip(inu,imu,jneigh,jatom))/y
            end do
          end do
        end if !end if y .gt. 1.0d-05)
      end do !ineig 
      imu = 0
      do issh = 1, nssh(in1)
        do mqn = 1, 2*lssh(issh,in1) + 1
          imu = imu + 1
          Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
        end do
        QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qout(issh,iatom)
      end do
    !Check whether there are negative charges and correct
    !If there's more than one shell whose charge is negative, more work is
    !needed, but that'd be quite pathological a situation...
    do issh = 1, nssh(in1)
      if( Qout(issh,iatom) .lt. 0 .and. nssh(in1) .gt. 1 ) then
        do jssh = 1,nssh(in1)
          if ( jssh .ne. issh ) then
            Qout(jssh,iatom) = Qout(jssh,iatom) + Qout(issh,iatom)/(nssh(in1)-1)
          end if !end if jssh .ne. issh 
        end do !end if jssh = 1,nssh(in1)
        Qout(issh,iatom) = 0.0d0              
      end if !end if  Qout(issh,iatom) .lt. 0
    end do !end do issh = 1, nssh(in1)
  end do !iatoms
end subroutine MULLIKEN_DIPOLE_CHARGES  


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

subroutine dipole_proyection()
  use M_system
  use M_constants
  use M_fdata, only: num_orb,nssh,lssh,Qneutral
  implicit none
  real, parameter ::  Debye = 0.208194
  real, parameter ::  klambda = 1.0
  integer iatom
  integer imu
  integer in1, in2
  integer ineigh
  integer inu
  integer issh
  integer jatom
  integer mbeta
  real      Qtot, Qtot1, Qtot2
  real, dimension(3) :: r1,r2,Rbc,u21
  real, dimension(3) :: rmedio, raux
  real      w_suma  
  real, dimension(3,3) :: bwrr, bwrr_inv, u_bwrr, ut_bwrr, v_bwrr, vt_bwrr, zero_bwrr
  real, dimension (natoms) :: c_k
  real, dimension (neigh_max) :: w_k
  real, dimension(3,natoms) :: intra_dip, res_dip
  real, dimension(3,1) :: intra_dip_aux, delta_ck
  integer n_bwrr, lda_bwrr, info, lda, i
  integer, dimension(3) :: ipiv
  real, dimension(3) :: lwork, s_bwrr, dummy
  do iatom = 1, natoms
    Q0_TOT(iatom) = 0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
    end do
  end do
  dip_x=0.0d0
  dip_y=0.0d0
  dip_z=0.0d0
  do iatom = 1, natoms
    Qtot=-Q0_TOT(iatom)
    in1 = imass(iatom)
    do issh = 1,nssh(in1)
      Qtot = Qtot+Qout(issh,iatom)
    end do
    dip_x = dip_x+Qtot*ratom(1,iatom)
    dip_y = dip_y+Qtot*ratom(2,iatom)
    dip_z = dip_z+Qtot*ratom(3,iatom)      
  enddo !end do iatom = 1,natoms
  dip_tot      = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )    
  dipQin_x    = dip_x
  dipQin_y    = dip_y
  dipQin_z    = dip_z
  dipQin_tot = dip_tot
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do issh = 1,nssh(in1)
     Qtot = Qtot+Qin(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    dip_x = dip_x-Q0_TOT(iatom)*r1(1)
    dip_y = dip_y-Q0_TOT(iatom)*r1(2)    
    dip_z = dip_z-Q0_TOT(iatom)*r1(3)
  end do !end do iatom = 1,natoms
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)    
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)      
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
  end do ! end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  dipTot_x    = dip_x
  dipTot_y    = dip_y
  dipTot_z    = dip_z
  dipTot_tot = dip_tot
  ! write the intraatomic dipole
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0      
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          if ((iatom .eq. jatom) .and. (imu .ne. inu)) then
            dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
            dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
            dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
          end if !end if
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
  end do ! end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  dipIntra_x    = dip_x
  dipIntra_y    = dip_y
  dipIntra_z    = dip_z
  dipIntra_tot = dip_tot
  !DIP RES = DIP_TOT - PROY
  dip_res_x = 0.0d0
  dip_res_y = 0.0d0
  dip_res_z = 0.0d0
  do iatom = 1, natoms
    dip_res_x = 0.0d0
    dip_res_y = 0.0d0
    dip_res_z = 0.0d0
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)    
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          if (iatom .ne. jatom) then !.or. imu .eq. inu) then
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1) + dipc(2,imu,inu,ineigh,iatom)*u21(2) + dipc(3,imu,inu,ineigh,iatom)*u21(3)
            dip_res_x = dip_res_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) - dip_proy*u21(1))
            dip_res_y = dip_res_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) - dip_proy*u21(2))
            dip_res_z = dip_res_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) - dip_proy*u21(3))
          end if !end if
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
    dip_res_tot = sqrt (dip_res_x**2 + dip_res_y**2 + dip_res_z**2 )
    dipRes_x = dipRes_x + dip_res_x
    dipRes_y = dipRes_y + dip_res_y
    dipRes_z = dipRes_z + dip_res_z
    res_dip(1,iatom) = dip_res_x
    res_dip(2,iatom) = dip_res_y
    res_dip(3,iatom) = dip_res_z
  end do ! end do iatom = 1,natoms
  dipRes_tot = sqrt (dipRes_x**2 + dipRes_y**2 + dipRes_z**2 ) 
  do iatom = 1, natoms
    dip_x = 0.0d0
    dip_y = 0.0d0
    dip_z = 0.0d0
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    Qtot=0.0d0
    Qtot1=0.0d0
    do issh = 1,nssh(in1)
      Qtot1 = Qtot1+Qin(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    Qtot=0.0d0
    Qtot2=0.0d0
    do issh = 1,nssh(in1)
      Qtot2 = Qtot2+Qout(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    jatom=iatom
    ineigh=neigh_self(iatom)
    in2=in1
    r2(:)=r1(:)
    Rbc(:)=(r1(:)+r2(:))/2.0d0
    do imu = 1,num_orb(in1)
      do inu = 1,num_orb(in2)
        if ( (imu .ne. inu)) then
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
        end if !end if
      end do !end do inu
    end do !end do imu
    intra_dip(1,iatom) = dip_x + res_dip(1,iatom)!/Debye
    intra_dip(2,iatom) = dip_y + res_dip(2,iatom)!/Debye
    intra_dip(3,iatom) = dip_z + res_dip(3,iatom)!/Debye
    dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  end do !end do iatom = 1,natoms
  dq_DP = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    rmedio = 0.0d0
    w_suma = 0.0d0
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      w_k(ineigh)=exp(-klambda*((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      rmedio = rmedio + w_k(ineigh)*r2
      w_suma = w_suma + w_k(ineigh)
    end do !end ineigh = 1,natoms
    rmedio = rmedio / w_suma
    bwrr = 0.0d0
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      do imu = 1,3 !xyz
        do inu = 1,3 !xyz
          bwrr(inu,imu) = bwrr(inu,imu) + w_k(ineigh)*r2(imu)*(r2(inu)-rmedio(inu))
        enddo !inu
      enddo !imu
    end do !end ineigh = 1,natoms
    !inversa de bwrr
    n_bwrr = 3
    lda_bwrr = n_bwrr
    lwork = n_bwrr
    bwrr_inv = bwrr
    call dgesvd('A', 'S', n_bwrr, n_bwrr, bwrr_inv, lda_bwrr,s_bwrr, u_bwrr, lda_bwrr, vt_bwrr, lda_bwrr, dummy, lwork, info)
    zero_bwrr = 0.00
    do i = 1,3
      if (abs(s_bwrr(i)) .gt. 0.000001) then
        zero_bwrr(i,i) = 1/s_bwrr(i)
      endif
    enddo
    v_bwrr = transpose(vt_bwrr)
    ut_bwrr = transpose(u_bwrr)
    bwrr_inv=matmul(v_bwrr,zero_bwrr)
    bwrr_inv=matmul(bwrr_inv,ut_bwrr)
    intra_dip_aux(1,1)=intra_dip(1,iatom)
    intra_dip_aux(2,1)=intra_dip(2,iatom)
    intra_dip_aux(3,1)=intra_dip(3,iatom)
    delta_ck = matmul(bwrr_inv,intra_dip_aux)
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      raux = r2(:) - rmedio(:)
      c_k(jatom) = w_k(ineigh) * (raux(1)*delta_ck(1,1)+raux(2)*delta_ck(2,1)+raux(3)*delta_ck(3,1))
      dq_DP(jatom) = dq_DP(jatom) + c_k(jatom)
    end do
  end do ! end do iatom = 1,natoms
  do iatom = 1, natoms
    Q0_TOT(iatom) = 0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
    end do
  end do
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0
  do iatom = 1, natoms
    Qtot=-Q0_TOT(iatom)+dq_DP(iatom)
    in1 = imass(iatom)
    do imu = 1,nssh(in1)
      Qtot = Qtot+Qout(imu,iatom)
    end do
    dip_x = dip_x+Qtot*ratom(1,iatom)
    dip_y = dip_y+Qtot*ratom(2,iatom)
    dip_z = dip_z+Qtot*ratom(3,iatom)
  end do !end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )    
  dipQout_x    = dip_x
  dipQout_y    = dip_y
  dipQout_z    = dip_z
  dipQout_tot = dip_tot
end subroutine dipole_proyection
  
real function get_min_parabola(x1,x2,x3,y1,y2,y3)
  real,intent(in) :: x1,x2,x3,y1,y2,y3
  real a,b,c
  b=(y2-y3)-(y1-y2)*(x2**2-x3**2)/(x1**2-x2**2)
  b=b/( (x2-x3)-(x1-x2)*(x2**2-x3**2)/(x1**2-x2**2) )
  a=(y1-y2)/(x1**2-x2**2)-b*(x1-x2)/(x1**2-x2**2)
  c=y1-a*x1**2-b*x1
  get_min_parabola=-b/(2*a)
  ! print*,a,'*x**2+',b,'*x+',c
end function get_min_parabola

