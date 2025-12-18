subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qin, Qout, &
  & rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, get_l_ofshell, &
  & get_orb_ofshell, get_issh_ofshell, g_h, g_xc, get_iatom_ofshell, ztot, &
  & g_h_shell, g_xc_shell
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer beta_iatom            
  integer imu, inu          
  integer in1, in2          
  integer issh, jssh
  integer ineigh ,jatom         
  integer mbeta
  real,dimension(nssh_tot,nssh_tot) :: A
  real,dimension(nssh_tot) :: c, SQ ! carga
  real :: Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real aux
  integer :: beta, alpha, alpha_iatom, ina, matom 
  integer :: alpha_iatom2, alpha2, nssh_tot2
  integer, dimension (:), allocatable :: mapindex
  integer, dimension (:), allocatable :: B
  integer, dimension (:,:), allocatable :: M
 
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
  call load_g_h_shell()
  call load_M()
  !Borrar es solo para hacer test----------
  !do issh=1,nsh_max
  !  fix_shell_charge(issh)=fix_shell_charge_aux(issh)
  !end do
  !lo usamos para H2O HsHsOsp
  !                   s s s p 
  fix_shell_charge = [1,0,1,0]
  print*,'fix_shell_charge',fix_shell_charge
  nssh_tot2 = count(fix_shell_charge == 0)
  allocate(mapindex(nssh_tot2))
  allocate(M(nssh_tot2+1,nssh_tot2+1))
  allocate(B(nssh_tot2+1))

 
  alpha_iatom2 = 0 
  
  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) alpha_iatom2=alpha_iatom2+1
    mapindex(alpha)=alpha_iatom2
  end do

  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) then
      do beta = 1, nssh_tot
        alpha_iatom  = get_iatom_ofshell(alpha)
        beta_iatom = get_iatom_ofshell(beta)
        matom = neigh_self(beta_iatom)
        if (fix_shell_charge(beta) .eq. 0) then
          print '(A,I0,A,I0,A,I0,A,I0,A)', 'M(', mapindex(alpha), ',', mapindex(beta), ') = f(', beta, ',', alpha, ')'
          M(mapindex(alpha),mapindex(beta)) = M(mapindex(alpha),mapindex(beta)) & 
                      & +  g_h_shell(beta , alpha ) &
                      & + g_xc_shell(beta , alpha ) &
                      & + g_xc_shell(alpha, beta  ) &
                      & - fshell_xc(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & - fshell_xc(alpha, beta  ,alpha_iatom, matom, beta_iatom) 
          do ineigh = 1, neighn(beta_iatom)
            mbeta = neigh_b(ineigh,beta_iatom)
            jatom = neigh_j(ineigh,beta_iatom)
            in2 = imass(jatom)
            in1 = imass(beta_iatom)
            do imu = 1, num_orb(in1)
              do inu = 1, num_orb(in2)
                aux = g_h(imu,inu,issh,alpha_iatom,ineigh,beta_iatom) + g_xc(imu,inu,issh,alpha_iatom,ineigh,beta_iatom)
                B(mapindex(alpha)) = B(mapindex(alpha)) + rho(imu,inu,ineigh,beta_iatom)*aux
              end do ! end do inu
            end do ! end do imu
            !B(mapindex(alpha)) = B(mapindex(alpha)) & 
            !           & + exc_shell_aa(alpha, alpha_iatom, neigh_beta_iatom, beta_iatom) &
            !           & - vxc_shell_aa(alpha, alpha_iatom, neigh_beta_iatom, beta_iatom)
          end do ! end do ineigh
        endif !fix_shell_charge(beta) = 0
        if (fix_shell_charge(beta) .eq. 1) then
          print '(A,I0,A,I0,A,I0,A,F8.4)', '*B*(', mapindex(alpha),') = f(', beta, ',', alpha, ')', Qin(get_issh_ofshell(beta),beta_iatom)
          B(mapindex(alpha)) = B(mapindex(alpha))- Qin(get_issh_ofshell(beta),beta_iatom)*( &
                      & +  g_h_shell(beta , alpha ) &
                      & + g_xc_shell(beta , alpha ) &
                      & + g_xc_shell(alpha, beta  ) &
                      & - fshell_xc(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & - fshell_xc(alpha, beta  ,alpha_iatom, matom, beta_iatom))!&
!                     & + exc_shell_aa(alpha)- vxc_shell_aa(alpha)  
          
        endif !fix_shell_charge(beta) = 1
      end do !beta
     end if !fix_shell_charge(alpha) = 0
    end do !alpha

    !SOLVE SYSTEM Mx = B.  x are the charges
    M(nssh_tot2+1,nssh_tot2+1) = 0
    B(nssh_tot2+1) = ztot

    do alpha = 1,nssh_tot
      if (fix_shell_charge(alpha) .eq. 0) then
        M(nssh_tot2+1,mapindex(alpha)) = 1
        M(mapindex(alpha),nssh_tot2+1) = 1
      else
        B(nssh_tot2+1) = B(nssh_tot2+1) - Qin(get_issh_ofshell(alpha),get_iatom_ofshell(alpha))
      endif
    end do
    !LWMAX = 100
    !call ssysv( 'U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, lwork, info )
    !call sgesv(nssh_tot,1,M,nssh_tot,ipiv,B,nssh_tot,info )
  stop    
  contains

    real function fshell_xc(alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom)
      use M_system, only: f_xc, get_orb_ofshell, get_l_ofshell, get_issh_ofshell
      implicit none
      integer, intent(in) :: alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom
      integer :: mu_min, mu_max, imu_local
      fshell_xc = 0.0
      mu_min = get_orb_ofshell(alpha)
      mu_max = mu_min + 2*get_l_ofshell(alpha)
      do imu_local = mu_min, mu_max
        fshell_xc = fshell_xc + f_xc(imu_local, get_issh_ofshell(beta), alpha_iatom, neigh_beta_iatom, beta_iatom)
      end do
      fshell_xc = fshell_xc /  (2*get_l_ofshell(alpha) + 1)
    end function fshell_xc

    real function exc_shell_aa(alpha, alpha_iatom, neigh_beta_iatom, beta_iatom)
      use M_system, only: exc_aa, get_orb_ofshell, get_l_ofshell, get_issh_ofshell
      implicit none
      integer, intent(in) :: alpha,alpha_iatom, neigh_beta_iatom, beta_iatom
      integer :: mu_min, mu_max, imu_local
      exc_shell_aa = 0.0
      mu_min = get_orb_ofshell(alpha)
      mu_max = mu_min + 2*get_l_ofshell(alpha)
      do imu_local = mu_min, mu_max
        !carlos exc_aa(imu)
        exc_shell_aa = exc_shell_aa + exc_aa(imu_local, alpha_iatom, neigh_beta_iatom, beta_iatom)
      end do
      exc_shell_aa = exc_shell_aa /  (2*get_l_ofshell(alpha) + 1)
    end function exc_shell_aa




    real function vxc_shell_aa(alpha, alpha_iatom, neigh_beta_iatom, beta_iatom)
      use M_system, only: vxc_aa, get_orb_ofshell, get_l_ofshell, get_issh_ofshell
      implicit none
      integer, intent(in) :: alpha,alpha_iatom, neigh_beta_iatom, beta_iatom
      integer :: mu_min, mu_max, imu_local
      vxc_shell_aa = 0.0
      mu_min = get_orb_ofshell(alpha)
      mu_max = mu_min + 2*get_l_ofshell(alpha)
      do imu_local = mu_min, mu_max
        vxc_shell_aa = vxc_shell_aa + vxc_aa(imu_local, alpha_iatom, neigh_beta_iatom, beta_iatom)
      end do
      vxc_shell_aa = vxc_shell_aa /  (2*get_l_ofshell(alpha) + 1)
    end function vxc_shell_aa

    subroutine load_M() !Mx=B
      use M_system, only: g_xc_shell,f_xc_shell,imass,vxc_aa_shell,exc_aa_shell,nssh_tot
      use M_fdata, only: gxc_1c,fxc_1c,exc_1c_0,vxc_1c_0, nssh
      implicit none
      integer :: iatom,count ,issh,kssh,alpha,beta
      g_xc_shell=0.0d0
      f_xc_shell=0.0d0
      gxc_1c=22.00
      count=0
      print*,'alpha,beta,count,issh,kssh,iatom  '
      do iatom=1,natoms
        do issh=1,nssh(imass(iatom))
          alpha = count + issh          
          do kssh=1, nssh(imass(iatom))
            beta = count + kssh
            print*,alpha,beta,count,issh,kssh,iatom  
            g_xc_shell(alpha,beta) = gxc_1c(imass(iatom),issh,issh,kssh)
            f_xc_shell(alpha,beta) = fxc_1c(imass(iatom),issh,kssh)
          end do
          exc_aa_shell(alpha) = exc_1c_0(imass(iatom),issh)
          vxc_aa_shell(alpha) = vxc_1c_0(imass(iatom),issh,issh)
        end do
        count = count + nssh(imass(iatom))
      end do

    do alpha = 1, nssh_tot
        write(*,'(4F10.2)') (g_xc_shell(alpha,beta), beta=1,nssh_tot)
    end do
     
    end subroutine 

    subroutine load_g_h_shell()
      use M_system, only: g_h_shell, natoms, get_shell_ofatom_imu, get_shell_ofatom_issh
      use M_fdata, only: nssh
      implicit none
      integer :: iatom, imu, alpha, matom, katom, kssh, beta
      g_h_shell = 0.0
      do iatom = 1, natoms
        do imu=1,  num_orb(imass(iatom))
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
    end subroutine  load_g_h_shell

end subroutine stationary_charges
