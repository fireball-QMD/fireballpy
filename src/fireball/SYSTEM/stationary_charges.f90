subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qin, Qout, &
  & rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, get_l_ofshell, &
  & get_orb_ofshell, get_issh_ofshell, g_h, g_xc, get_iatom_ofshell, ztot, &
  & g_h_shell, g_xc_shell,f_xc_shell,exc_aa_shell,vxc_aa_shell
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer imu, inu          
  integer in1, in2          
  integer issh, jssh, alpha_issh, beta_iatom
  integer ineigh ,jatom         
  integer i,j
  real(double),dimension(nssh_tot,nssh_tot) :: A
  real(double),dimension(nssh_tot) :: c, SQ ! carga
  real(double) :: Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real(double) aux
  integer :: beta, alpha, alpha_iatom, ina, matom, iatom 
  integer :: alpha_iatom2, alpha2, nssh_tot2 
  integer, dimension (:), allocatable :: mapindex, ipiv
  integer :: info, lwork 
  real(double), dimension (:), allocatable :: B, work
  real(double), dimension (:,:), allocatable :: M
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
  call load_M()
  !Borrar es solo para hacer test----------
  !do issh=1,nsh_max
  !  fix_shell_charge(issh)=fix_shell_charge_aux(issh)
  !end do
  !lo usamos para H2O HsHsOsp
  !                   s s s p 
  !fix_shell_charge = [1,0,1,0]
  !                   CCC
  fix_shell_charge = [0,0,0,0,0,0]

  print*,'fix_shell_charge',fix_shell_charge
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
  B(nssh_tot2+1) = ztot
 
  alpha_iatom2 = 0 
  
  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) alpha_iatom2=alpha_iatom2+1
    mapindex(alpha)=alpha_iatom2
  end do

  print*, '-----Qin----'
  print*, Qin

  print*, '-----g_h----'
  do i = 1, nssh_tot
      write(*,'(100(1x,ES14.6))') (g_h_shell(i,j), j=1,nssh_tot)
  end do
  print*, '-----g_xc----'
  do i = 1, nssh_tot
      write(*,'(100(1x,ES14.6))') (g_xc_shell(i,j), j=1,nssh_tot)
  end do

  print*, '-----f_xc----'
  do i = 1, nssh_tot
      write(*,'(100(1x,ES14.6))') (f_xc_shell(i,j), j=1,nssh_tot)
  end do



  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) /= 0) cycle
    do beta = 1, nssh_tot
      if (fix_shell_charge(beta) .eq. 0) then
        M(mapindex(alpha),mapindex(beta)) = M(mapindex(alpha),mapindex(beta)) & 
        & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
        & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  )
      else
        B(mapindex(alpha)) = B(mapindex(alpha))- Qin(get_issh_ofshell(beta),beta_iatom)*( &
        & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
        & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  ))!&
      endif !fix_shell_charge(beta) = 1
    end do !beta
    alpha_iatom = get_iatom_ofshell(alpha)
    alpha_issh = get_issh_ofshell(alpha)
    do iatom = 1, natoms
      in1 = imass(iatom) 
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        do imu = 1, num_orb(in1)
          do inu = 1, num_orb(in2)
            aux = 0.0d0*g_h(imu,inu,alpha_issh,alpha_iatom,ineigh,iatom) + g_xc(imu,inu,alpha_issh,alpha_iatom,ineigh,iatom)
            B(mapindex(alpha)) = B(mapindex(alpha)) +  rho(imu,inu,ineigh,iatom)*aux
          end do ! inu
        end do ! imu
      end do ! ineigh
    end do ! iatom
    B(mapindex(alpha)) = B(mapindex(alpha)) + exc_aa_shell(alpha) - vxc_aa_shell(alpha)
  end do !alpha


!    do alpha = 1,nssh_tot
!      if (fix_shell_charge(alpha) .eq. 0) then
!        M(nssh_tot2+1,mapindex(alpha)) = 1
!        M(mapindex(alpha),nssh_tot2+1) = 1
!      else
!        B(nssh_tot2+1) = B(nssh_tot2+1) - Qin(get_issh_ofshell(alpha),get_iatom_ofshell(alpha))
!      endif
!    end do

    print *, "===== MATRIZ M ====="
    do i = 1, nssh_tot2+1
        write(*,'(100(1x,ES14.6))') (M(i,j), j=1,nssh_tot2+1)
    end do

    print *, "===== VECTOR B ====="
    do i = 1, nssh_tot2+1
        write(*,'(1x,ES14.6)') B(i)
    end do


    allocate(ipiv(nssh_tot))
    allocate(work(1))
    call dsysv('U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, -1, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsysv('U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, lwork, info)
    deallocate(work)
    deallocate(ipiv)
    if (info /= 0) then
       print *, "Error in dsysv, info =", info
    end if
   
    do alpha = 1, nssh_tot
       issh=get_issh_ofshell(alpha)
       iatom=get_iatom_ofshell(alpha)
       print*,'stationary_charges',issh,iatom,mapindex(alpha),fix_shell_charge(alpha)
       if (fix_shell_charge(alpha) .eq. 0) then
          Qout(issh,iatom) = B(mapindex(alpha))
       else
          Qout(issh,iatom) = Qin(issh,iatom)
       end if
    end do

    print *, "===== VECTOR B out ====="
    do i = 1, nssh_tot2+1
        write(*,'(1x,ES14.6)') B(i)
    end do
    stop
  contains

    subroutine load_M() !Mx=B
      use M_system, only: g_xc_shell,f_xc_shell,imass,vxc_aa_shell,exc_aa_shell,nssh_tot, &
      & g_h_shell, natoms, get_shell_ofatom_imu, get_shell_ofatom_issh
      use M_fdata, only: gxc_1c,fxc_1c,exc_1c_0,vxc_1c_0, nssh, Qneutral
      implicit none
      integer :: iatom,count ,issh,kssh,alpha,beta,imu, matom, katom 
      g_h_shell = 0.0d0
      g_xc_shell=0.0d0
      vxc_aa_shell = 0.0d0
      f_xc_shell=0.0d0
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
              & g_h(imu,imu, kssh, katom, matom, iatom) / (2*get_l_ofshell(alpha) + 1)
            end do
          end do
        end do
      end do
    end subroutine  load_M

end subroutine stationary_charges
