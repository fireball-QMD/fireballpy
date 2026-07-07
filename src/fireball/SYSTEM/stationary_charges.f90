subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qin, Qout, &
  & rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, get_l_ofshell, &
  & get_orb_ofshell, get_issh_ofshell, g_h, g_xc, get_iatom_ofshell, ztot, &
  & g_h_shell, g_xc_shell,f_xc_shell,exc_aa_shell,vxc_aa_shell, qstate, symbol, &
  & get_shell_ofatom_issh, Kscf, ratom
  use M_fdata, only: num_orb,nssh,lssh,Qneutral
  implicit none
  integer imu, inu          
  integer in1, in2          
  integer issh, jssh, beta_iatom
  integer ineigh ,jatom         
  integer i,j
  real(double),dimension(nssh_tot,nssh_tot) :: A
  real(double),dimension(nssh_tot) :: c, SQ ! carga
  ! correcciones GSN del sistema completo linealizado (igsn=4), calculadas en load_M
  real(double),dimension(nssh_tot,nssh_tot) :: gsnM
  real(double),dimension(nssh_tot) :: gsnB
  real(double) :: Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real(double) aux, ztot_aux
  integer :: beta, alpha, ina, matom, iatom
  integer :: alpha_iatom2, alpha2, nssh_tot2 
  integer, dimension (:), allocatable :: mapindex, ipiv
  integer :: info, lwork
  real(double), dimension (:), allocatable :: B, work
  real(double), dimension (:,:), allocatable :: M
  ! isolver: 0 = dgesv directo; 1 = dgelsd (truncamiento SVD; OJO: saltos de rango entre
  ! geometrias -> PES discontinua, solo para diagnostico); 2 = Tikhonov suave en dQ=Q-Q0
  ! con ligadura exacta (recomendado con shells casi degeneradas: 2+ atomos O, L repetido).
  integer, parameter :: isolver = 0
  real(double), parameter :: rcond_sc = 1.0d-3   ! umbral relativo isolver=1
  real(double), parameter :: ridge_sc = 3.0d-2   ! lambda de Tikhonov isolver=2 (unidades de M ~ eV/e)
  integer, parameter :: ifix_dshell = 2          ! 1 = shells d fijas a carga neutra (ver abajo)
  integer :: rank_sc
  integer, dimension (:), allocatable :: iwork_lsd
  real(double), dimension (:), allocatable :: Q0v, sing, bp, xv
  real(double), dimension (:,:), allocatable :: Gm
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
  call load_M()
  !todos libres
  do issh=1,nssh_tot
    fix_shell_charge(issh)=0
  end do
  ! ifix_dshell=1: fijar las shells d (polarizacion vacia, Q0=0) a su carga neutra.
  ! ifix_dshell=2: fijar TODAS las shells con Qneutral=0 (d y excitadas s*/p* de base doble).
  ! Motivo: las shells vacias libres crean direcciones casi nulas del sistema estacionario
  ! (dimero de agua: Q_d=-1.0, gap 10 mA / 0.55 eV/A; base doble: SCF no converge, Q basura);
  ! fijandolas el sistema es variacional (dimero 0.03 mA / 0.004 eV/A).
  if (ifix_dshell .eq. 1) then
    do issh=1,nssh_tot
      if (get_l_ofshell(issh) .eq. 2) fix_shell_charge(issh)=1
    end do
  else if (ifix_dshell .eq. 2) then
    do issh=1,nssh_tot
      if (Qneutral(get_issh_ofshell(issh), imass(get_iatom_ofshell(issh))) .lt. 1.0d-10) fix_shell_charge(issh)=1
    end do
  end if
  !fix_shell_charge(1)=1
  !lo usamos para H2O HsHsOsp
  !                   s s s p 
  !fix_shell_charge = [1,0,1,0]
  !                   CC
  !fix_shell_charge = [1,0,1,1,0,1]
  !fix_shell_charge = [0,1,1,0,1,1] 
  !fix_shell_charge = [1,0,1,1,0,1,0,1] 
  !fix_shell_charge = [1,0,1,0]
  ztot_aux=0.0d0
  !ztot_aux=qstate

  do alpha=1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) ztot_aux = ztot_aux + Qin(get_issh_ofshell(alpha),get_iatom_ofshell(alpha))
  end do

  if (Kscf .eq. 1) then
    print*, "========  positions ========",'Kscf = ',Kscf
    do iatom = 1, natoms
       in1 = imass(iatom)
       write (*,'(2X,A2,3F12.6)') symbol(iatom), ratom(:,iatom)  - (/ 3.141593, 0.367879, 1.414214 /)
     enddo
   end if
   print*, "======== Qin  CHARGES ======== Kscf = " ,Kscf
   do iatom = 1, natoms
     in1 = imass(iatom)
     aux = sum(Qin(1:nssh(in1), iatom))
     write(*,'(2X,A2,1X," | ",F12.6," |",100(1X,F12.6,1X,"(",I0,")"," |"))') &
     symbol(iatom), aux, &
     (Qin(issh,iatom), &
     fix_shell_charge(get_shell_ofatom_issh(iatom,issh)), &
     issh = 1, nssh(in1))
   end do
   print*, "=============================="


  print*,'ztot_fix =',ztot-ztot_aux
  nssh_tot2 = count(fix_shell_charge == 0)
  allocate(mapindex(nssh_tot2))
  allocate(M(nssh_tot2+1,nssh_tot2+1))
  allocate(B(nssh_tot2+1))
  M=0.0d0
  B=0.0d0
  !SOLVE SYSTEM Mx = B.  x are the charges
  M(:,nssh_tot2+1) = 1.0
  M(nssh_tot2+1,:) = 1.0
  M(nssh_tot2+1,nssh_tot2+1)= 0.0
  B(nssh_tot2+1) = ztot_aux
 
  alpha_iatom2 = 0 
  
  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) alpha_iatom2=alpha_iatom2+1
    mapindex(alpha)=alpha_iatom2
  end do


  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) /= 0) cycle
    do beta = 1, nssh_tot
      if (fix_shell_charge(beta) .eq. 0) then
        M(mapindex(alpha),mapindex(beta)) = M(mapindex(alpha),mapindex(beta)) &
        & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
        & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  ) + gsnM(alpha,beta)
      else
        ! (bug corregido: beta_iatom no estaba inicializado)
        B(mapindex(alpha)) = B(mapindex(alpha))- Qin(get_issh_ofshell(beta),get_iatom_ofshell(beta))*( &
        & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
        & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  ) + gsnM(alpha,beta))!&
      endif !fix_shell_charge(beta) = 1
    end do !beta
    do iatom = 1, natoms
      in1 = imass(iatom) 
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        do imu = 1, num_orb(in1)
          do inu = 1, num_orb(in2)
            aux = g_h(alpha,imu,inu,ineigh,iatom) + g_xc(alpha,imu,inu,ineigh,iatom)
            B(mapindex(alpha)) = B(mapindex(alpha)) +  rho(imu,inu,ineigh,iatom)*aux
          end do ! inu
        end do ! imu
      end do ! ineigh
    end do ! iatom
    B(mapindex(alpha)) = B(mapindex(alpha)) + exc_aa_shell(alpha) - vxc_aa_shell(alpha) + gsnB(alpha)
  end do !alpha


  call print_matrix_as_numpy("M", M, nssh_tot2, nssh_tot2)
  call print_vector_as_numpy("B", B, nssh_tot2)
    
    ! simetrizo M = (M + M^T)/2
    !do i = 1, nssh_tot2
    !  do j = i+1, nssh_tot2
    !    aux = (M(i,j) + M(j,i)) / 2.0d0
    !    M(i,j) = aux
    !    M(j,i) = aux
    !  end do
    !end do

  if (isolver .eq. 0) then
    ! solve directo (sin regularizar)
    allocate(ipiv(nssh_tot2+ 1))
    call dgesv(nssh_tot2+ 1, 1, M, nssh_tot2+ 1, ipiv, B, nssh_tot2+ 1, info)
    deallocate(ipiv)
    if (info /= 0) print *, "Error in dgesv, info =", info
  else if (isolver .eq. 1) then
    ! solve regularizado SVD (dgelsd) en dQ = Q - Q0: las direcciones casi nulas
    ! (valor singular < rcond_sc * sigma_max) quedan con dQ = 0 (cargas neutras)
    allocate(Q0v(nssh_tot2+1))
    Q0v = 0.0d0
    do alpha = 1, nssh_tot
      if (fix_shell_charge(alpha) .eq. 0) then
        Q0v(mapindex(alpha)) = Qneutral(get_issh_ofshell(alpha), imass(get_iatom_ofshell(alpha)))
      end if
    end do
    B = B - matmul(M, Q0v)
    allocate(sing(nssh_tot2+1))
    allocate(work(1), iwork_lsd(1))
    call dgelsd(nssh_tot2+1, nssh_tot2+1, 1, M, nssh_tot2+1, B, nssh_tot2+1, sing, rcond_sc, rank_sc, work, -1, iwork_lsd, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    deallocate(iwork_lsd)
    allocate(iwork_lsd(max(1, 30*(nssh_tot2+1))))
    call dgelsd(nssh_tot2+1, nssh_tot2+1, 1, M, nssh_tot2+1, B, nssh_tot2+1, sing, rcond_sc, rank_sc, work, lwork, iwork_lsd, info)
    deallocate(work, iwork_lsd)
    if (info /= 0) print *, "Error in dgelsd, info =", info
    if (Kscf .eq. 1) then
      write (*,'(A,I4,A,I4,A,ES10.3,A,ES10.3)') ' dgelsd: rank=', rank_sc, ' /', nssh_tot2+1, &
        & '  smax=', sing(1), '  smin=', sing(nssh_tot2+1)
    end if
    deallocate(sing)
    B = B + Q0v   ! Q = Q0 + dQ  (la componente n+1 es mu, con Q0v(n+1)=0)
    deallocate(Q0v)
  end if

  if (isolver .eq. 2) then
    ! Tikhonov suave en dQ = Q - Q0, ligadura exacta, mu sin penalizar:
    !   min ||M_fis dQ + mu*1 - b'||^2 + ridge_sc^2 ||dQ||^2   s.a.  sum(dQ) = qlibre
    ! (suave con la geometria: sin saltos de rango como el truncamiento SVD)
    allocate(Q0v(nssh_tot2+1))
    Q0v = 0.0d0
    do alpha = 1, nssh_tot
      if (fix_shell_charge(alpha) .eq. 0) then
        Q0v(mapindex(alpha)) = Qneutral(get_issh_ofshell(alpha), imass(get_iatom_ofshell(alpha)))
      end if
    end do
    allocate(bp(nssh_tot2), Gm(nssh_tot2+2,nssh_tot2+2), xv(nssh_tot2+2))
    bp = B(1:nssh_tot2) - matmul(M(1:nssh_tot2,1:nssh_tot2), Q0v(1:nssh_tot2))
    Gm = 0.0d0
    ! A1 = [M_fis | 1] (n x n+1); G = A1^T A1 + ridge^2 (solo bloque dQ)
    Gm(1:nssh_tot2+1,1:nssh_tot2+1) = matmul(transpose(M(1:nssh_tot2,1:nssh_tot2+1)), M(1:nssh_tot2,1:nssh_tot2+1))
    do i = 1, nssh_tot2
      Gm(i,i) = Gm(i,i) + ridge_sc*ridge_sc
    end do
    Gm(1:nssh_tot2,nssh_tot2+2) = 1.0d0
    Gm(nssh_tot2+2,1:nssh_tot2) = 1.0d0
    xv(1:nssh_tot2+1) = matmul(transpose(M(1:nssh_tot2,1:nssh_tot2+1)), bp)
    xv(nssh_tot2+2) = B(nssh_tot2+1) - sum(Q0v(1:nssh_tot2))
    allocate(ipiv(nssh_tot2+2))
    call dgesv(nssh_tot2+2, 1, Gm, nssh_tot2+2, ipiv, xv, nssh_tot2+2, info)
    deallocate(ipiv)
    if (info /= 0) print *, "Error in dgesv (tikhonov), info =", info
    B(1:nssh_tot2) = Q0v(1:nssh_tot2) + xv(1:nssh_tot2)
    B(nssh_tot2+1) = xv(nssh_tot2+1)   ! mu
    deallocate(bp, Gm, xv, Q0v)
  end if

!    allocate(ipiv(nssh_tot + 1))
!    allocate(work(1))
!    call dsysv('U', nssh_tot + 1, 1, M, nssh_tot + 1, ipiv, B, nssh_tot + 1, work, -1, info)
!    lwork = int(work(1))
!    deallocate(work)
!    allocate(work(lwork))
!    call dsysv('U', nssh_tot + 1, 1, M, nssh_tot + 1, ipiv, B, nssh_tot + 1, work, lwork, info)
!    deallocate(work)
!    deallocate(ipiv)
!    if (info /= 0) then
!       print *, "Error in dsysv, info =", info
!    end if
!
    print *, "===== VECTOR B out ====="
    do i = 1, nssh_tot2+1
        write(*,'(1x,F14.6)') B(i)
    end do
    aux = 0.0d0
    do i = 1, nssh_tot2
      aux = aux + B(i)
    end do
    print*,'sum B = ',aux
    
    do alpha = 1, nssh_tot
       issh=get_issh_ofshell(alpha)
       iatom=get_iatom_ofshell(alpha)
       if (fix_shell_charge(alpha) .eq. 0) then
          Qout(issh,iatom) = B(mapindex(alpha))
       else
          Qout(issh,iatom) = Qin(issh,iatom)
       end if
    end do

   print*, "======== Qout CHARGES ======== Kscf = " ,Kscf
   do iatom = 1, natoms
     in1 = imass(iatom)
     aux = sum(Qout(1:nssh(in1), iatom))
     write(*,'(2X,A2,1X," | ",F12.6," |",100(1X,F12.6,1X,"(",I0,")"," |"))') &
     symbol(iatom), aux, &
     (Qout(issh,iatom), &
     fix_shell_charge(get_shell_ofatom_issh(iatom,issh)), &
     issh = 1, nssh(in1))
   end do
   print*, " =============================="

!   stop
  contains

    subroutine load_M() !Mx=B
      use M_system, only: g_xc_shell,f_xc_shell,imass,vxc_aa_shell,exc_aa_shell,nssh_tot, &
      & g_h_shell, natoms, get_shell_ofatom_imu, get_shell_ofatom_issh, den_sh, den_or, &
      & arho_on, arhoi_on, rho_on, rhoi_on, arho_off, arhoij_off, rho_off, rhoij_off, &
      & s_mat, sm_mat, xc_overtol, orb2shell, neighn, neigh_j, neigh_b
      use M_fdata, only: gxc_1c,fxc_1c,exc_1c_0,vxc_1c_0, nssh, Qneutral
      implicit none
      integer :: iatom,count ,issh,kssh,alpha,beta,imu, matom, katom
      integer :: kmatom, kalpha
      integer :: ineigh, jatom, in1, in2, jssh, inu
      logical :: onsite
      ! igsn: 0=sin GSN; 1=ec.117 (simplificado, congelado en Q^R); 2=ec.117+Delta_alpha (ec.119);
      !       3=ec.117-Delta_alpha (test de signo); 4=sistema COMPLETO linealizado (III.D, opcion (b):
      !       coeficientes vivos en Qin -> al converger resuelve la estacionariedad exacta de la energia viva)
      ! Medido en scan fino H2O (fdata_superspline), gap E-min vs F=0 / max|Fz+dE/dz|:
      ! igsn=0: 1.21 mA / 0.055 eV/A ; igsn=1: 3.22 mA / 0.35 ; igsn=2: 1.83 mA / 0.43 ; igsn=3: 5.16 mA / 0.30
      ! igsn=4: 0.22 mA / 0.012 eV/A  (mulliken_dipole como referencia: 3.02 mA / 0.15)
      integer, parameter :: igsn = 4
      real(double) :: nR_aa, naR_aa, qk, exc_R, muxc_R, exc_aR, muxc_aR
      real(double) :: dexc_dum, d2exc_dum, dmuxc_dum, d2muxc_dum
      real(double) :: Mba, nR_bb, naR_bb, qg, dexc_b, dmuxc_b, exc_dum, muxc_dum, delta_a
      real(double) :: dexc_b2, dmuxc_b2, Ssh, abar, abarL, dmu, d2mu, dmuL, d2muL
      real(double) :: dn, dnL, Smunu, dv, pref, Mv, Nv
      g_h_shell = 0.0d0
      g_xc_shell=0.0d0
      vxc_aa_shell = 0.0d0
      f_xc_shell=0.0d0
      gsnM = 0.0d0
      gsnB = 0.0d0
      count=0
      do iatom=1,natoms
        do issh=1,nssh(imass(iatom))
          alpha = count + issh          
          exc_aa_shell(alpha) = exc_1c_0(issh,imass(iatom))
          vxc_aa_shell(alpha) = vxc_1c_0(issh,issh,imass(iatom))
          do kssh=1, nssh(imass(iatom))
            beta = count + kssh
            f_xc_shell(alpha,beta) = fxc_1c(issh,kssh,imass(iatom))
            g_xc_shell(alpha,beta) = gxc_1c(issh,issh,kssh,imass(iatom))
            exc_aa_shell(alpha) = exc_aa_shell(alpha) - Qneutral(kssh, imass(iatom))*f_xc_shell(alpha,beta)
            vxc_aa_shell(alpha) = vxc_aa_shell(alpha) - Qneutral(kssh, imass(iatom))*g_xc_shell(alpha,beta)
          end do
        end do
        count = count + nssh(imass(iatom))
      end do
      do iatom = 1, natoms
        do imu=1,num_orb(imass(iatom))
          alpha = get_shell_ofatom_imu(iatom,imu)
          matom=neigh_self(iatom)
          do katom = 1, natoms
            do kssh = 1, nssh(imass(katom))
              beta = get_shell_ofatom_issh(katom,kssh)
              g_h_shell(alpha,beta) = g_h_shell(alpha,beta) +&
              & g_h(beta,imu,imu, matom, iatom) / (2*get_l_ofshell(alpha) + 1)
            end do
          end do
        end do
      end do
      ! GSN terms of the simplified stationary system (eq. 117, XC_basic doc):
      ! eps_xc(n^R_aa) - eps_xc(n^aR_aa) - V_xc(n^R_aa) + V_xc(n^aR_aa),
      ! frozen at the reference (neutral) charges. n^R_aa = sum_gamma Q0_gamma M^gamma_aa
      ! via den_sh; n^aR_aa restricts gamma to shells of the same atom.
      ! igsn = 1..3: version simplificada ec. 117 (congelada en Q^R), con/sin Delta_alpha (ec. 119)
      if (igsn .ge. 1 .and. igsn .le. 3) then
      do alpha = 1, nssh_tot
        iatom = get_iatom_ofshell(alpha)
        issh = get_issh_ofshell(alpha)
        matom = neigh_self(iatom)
        nR_aa = 0.0d0
        naR_aa = 0.0d0
        do beta = 1, nssh_tot
          katom = get_iatom_ofshell(beta)
          qk = Qneutral(get_issh_ofshell(beta), imass(katom))
          nR_aa = nR_aa + qk*den_sh(beta,issh,issh,matom,iatom)
          if (katom .eq. iatom) naR_aa = naR_aa + qk*den_sh(beta,issh,issh,matom,iatom)
        end do
        call cepal (nR_aa, exc_R, muxc_R, dexc_dum, d2exc_dum, dmuxc_dum, d2muxc_dum)
        call cepal (naR_aa, exc_aR, muxc_aR, dexc_dum, d2exc_dum, dmuxc_dum, d2muxc_dum)
        exc_aa_shell(alpha) = exc_aa_shell(alpha) + exc_R - exc_aR
        vxc_aa_shell(alpha) = vxc_aa_shell(alpha) + muxc_R - muxc_aR
        ! Delta_alpha (eq. 119): GSN response surviving at Q=Q^R,
        ! Delta = sum_beta M^beta_aa (e'-V')(n^R_bb) Q^R_beta
        !       - sum_{beta in a} M^beta_aa (e'-V')(n^{a,R}_bb) Q^R_beta
        if (igsn .ge. 2) then
          delta_a = 0.0d0
          do beta = 1, nssh_tot
            katom = get_iatom_ofshell(beta)
            kssh = get_issh_ofshell(beta)
            qk = Qneutral(kssh, imass(katom))
            Mba = den_sh(beta,issh,issh,matom,iatom)
            if (abs(Mba*qk) .lt. 1.0d-12) cycle
            kmatom = neigh_self(katom)
            nR_bb = 0.0d0
            naR_bb = 0.0d0
            do kalpha = 1, nssh_tot
              qg = Qneutral(get_issh_ofshell(kalpha), imass(get_iatom_ofshell(kalpha)))
              nR_bb = nR_bb + qg*den_sh(kalpha,kssh,kssh,kmatom,katom)
              if (get_iatom_ofshell(kalpha) .eq. katom) naR_bb = naR_bb + qg*den_sh(kalpha,kssh,kssh,kmatom,katom)
            end do
            call cepal (nR_bb, exc_dum, muxc_dum, dexc_b, d2exc_dum, dmuxc_b, d2muxc_dum)
            delta_a = delta_a + Mba*(dexc_b - dmuxc_b)*qk
            if (katom .eq. iatom) then
              call cepal (naR_bb, exc_dum, muxc_dum, dexc_b, d2exc_dum, dmuxc_b, d2muxc_dum)
              delta_a = delta_a - Mba*(dexc_b - dmuxc_b)*qk
            end if
          end do
          if (igsn .eq. 2) exc_aa_shell(alpha) = exc_aa_shell(alpha) + delta_a
          if (igsn .eq. 3) exc_aa_shell(alpha) = exc_aa_shell(alpha) - delta_a
        end if
      end do
      end if

      ! igsn = 4: sistema COMPLETO linealizado (seccion III.D del doc, opcion (b): los
      ! coeficientes epsilon'/V' se evaluan en las densidades vivas construidas con Qin;
      ! al converger el SCF se satisface la estacionariedad exacta de la energia viva).
      if (igsn .eq. 4) then
        ! (i) constantes GSN vivas de (eps_aa - v_aa), ecs. 101-102 evaluadas en Qin:
        !     eps_xc(n_aa) - eps_xc(n^a_aa) - V_xc(n_aa) + V_xc(n^a_aa)
        do alpha = 1, nssh_tot
          iatom = get_iatom_ofshell(alpha)
          issh = get_issh_ofshell(alpha)
          call cepal (arho_on(issh,issh,iatom), exc_R, muxc_R, dexc_dum, d2exc_dum, dmuxc_dum, d2muxc_dum)
          call cepal (arhoi_on(issh,issh,iatom), exc_aR, muxc_aR, dexc_dum, d2exc_dum, dmuxc_dum, d2muxc_dum)
          exc_aa_shell(alpha) = exc_aa_shell(alpha) + exc_R - exc_aR
          vxc_aa_shell(alpha) = vxc_aa_shell(alpha) + muxc_R - muxc_aR
        end do
        ! (ii) respuesta GSN del doble conteo en la MATRIZ (ecs. 103-106):
        !     A(alpha,beta) += -(e'-V')(n_bb) dn_bb/dQ_alpha [+ (e'-V')(n^b_bb) dn^b_bb/dQ_alpha si alpha en b]
        !     con dn_bb/dQ_alpha = M^alpha_bb / S_bb (den_sh / solape esferico on-site)
        do beta = 1, nssh_tot
          katom = get_iatom_ofshell(beta)
          kssh = get_issh_ofshell(beta)
          kmatom = neigh_self(katom)
          Ssh = sm_mat(kssh,kssh,kmatom,katom)
          if (abs(Ssh) .lt. xc_overtol) Ssh = sign(xc_overtol, Ssh)
          call cepal (arho_on(kssh,kssh,katom), exc_dum, muxc_dum, dexc_b, d2exc_dum, dmuxc_b, d2muxc_dum)
          call cepal (arhoi_on(kssh,kssh,katom), exc_dum, muxc_dum, dexc_b2, d2exc_dum, dmuxc_b2, d2muxc_dum)
          do alpha = 1, nssh_tot
            Mba = den_sh(alpha,kssh,kssh,kmatom,katom)/Ssh
            if (abs(Mba) .lt. 1.0d-14) cycle
            gsnM(alpha,beta) = gsnM(alpha,beta) - (dexc_b - dmuxc_b)*Mba
            if (get_iatom_ofshell(alpha) .eq. katom) gsnM(alpha,beta) = gsnM(alpha,beta) + (dexc_b2 - dmuxc_b2)*Mba
          end do
        end do
        ! (iii) respuesta de banda (ec. 91): gsnB(alpha) = sum_{pares,munu} P_numu dv^olsxc_munu/dQ_alpha
        !     dv/dQ = V''(nbar)(M^a/S)(n_munu - S_munu nbar) + V'(nbar) N^a
        !           - [V''(nbarL)(M^a/S)(nL_munu - S_munu nbarL) + V'(nbarL) N^a]  (solo alpha en a,b)
        do iatom = 1, natoms
          in1 = imass(iatom)
          do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            onsite = (iatom .eq. jatom .and. neigh_b(ineigh,iatom) .eq. 0)
            do issh = 1, nssh(in1)
              do jssh = 1, nssh(in2)
                if (onsite) then
                  abar  = arho_on(issh,jssh,iatom)
                  abarL = arhoi_on(issh,jssh,iatom)
                else
                  abar  = arho_off(issh,jssh,ineigh,iatom)
                  abarL = arhoij_off(issh,jssh,ineigh,iatom)
                end if
                Ssh = sm_mat(issh,jssh,ineigh,iatom)
                if (abs(Ssh) .lt. xc_overtol) Ssh = sign(xc_overtol, Ssh)
                call cepal (abar,  exc_dum, muxc_dum, dexc_dum, d2exc_dum, dmu,  d2mu)
                call cepal (abarL, exc_dum, muxc_dum, dexc_dum, d2exc_dum, dmuL, d2muL)
                do imu = 1, num_orb(in1)
                  if (orb2shell(imu,in1) .ne. issh) cycle
                  do inu = 1, num_orb(in2)
                    if (orb2shell(inu,in2) .ne. jssh) cycle
                    if (onsite) then
                      dn  = rho_on(imu,inu,iatom)
                      dnL = rhoi_on(imu,inu,iatom)
                      Smunu = 0.0d0
                      if (imu .eq. inu) Smunu = 1.0d0
                    else
                      dn  = rho_off(imu,inu,ineigh,iatom)
                      dnL = rhoij_off(imu,inu,ineigh,iatom)
                      Smunu = s_mat(imu,inu,ineigh,iatom)
                    end if
                    pref = rho(imu,inu,ineigh,iatom)
                    if (abs(pref) .lt. 1.0d-14) cycle
                    do alpha = 1, nssh_tot
                      Mv = den_sh(alpha,issh,jssh,ineigh,iatom)/Ssh
                      Nv = den_or(alpha,imu,inu,ineigh,iatom)
                      if (abs(Mv) .lt. 1.0d-14 .and. abs(Nv) .lt. 1.0d-14) cycle
                      dv = d2mu*Mv*(dn - Smunu*abar) + dmu*Nv
                      katom = get_iatom_ofshell(alpha)
                      if (katom .eq. iatom .or. katom .eq. jatom) then
                        dv = dv - d2muL*Mv*(dnL - Smunu*abarL) - dmuL*Nv
                      end if
                      gsnB(alpha) = gsnB(alpha) + pref*dv
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
    end subroutine  load_M
  subroutine print_matrix_as_numpy(name, M, nrows, ncols)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: nrows, ncols
    real(8), intent(in) :: M(nrows+1, ncols+1)
    integer :: i, j

    write(*,'(A)', advance='no') trim(name) // " = np.array(["
    
    do i = 1, nrows
        write(*,'(A)', advance='no') "    ["
        
        do j = 1, ncols
            if (j < ncols) then
                write(*,'(F12.6,A)', advance='no') M(i,j), ", "
            else
                write(*,'(F12.6)', advance='no') M(i,j)
            end if
        end do
        
        if (i < nrows) then
            print *, "],"
        else
            print *, "]"
        end if
    end do

    print *, "])"
    
  end subroutine print_matrix_as_numpy
  subroutine print_vector_as_numpy(name, v, n)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(8), intent(in) :: v(n+1)
    integer :: i
    write(*,'(A)', advance='no') trim(name) // " = np.array(["
    do i = 1, n
        if (i < n) then
            write(*,'(F12.6,A)', advance='no') v(i), ", "
        else
            write(*,'(F12.6)', advance='no') v(i)
        end if
    end do
    print *, "])"
  end subroutine print_vector_as_numpy
end subroutine stationary_charges
