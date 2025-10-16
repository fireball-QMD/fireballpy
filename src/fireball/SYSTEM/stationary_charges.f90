subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qin, Qout, &
  & rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, get_l_ofshell, &
  & get_orb_ofshell, get_issh_ofshell, g_h, g_xc, get_iatom_ofshell, ztot
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
  real :: diff_err,Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real aux
  integer :: beta, alpha, alpha_iatom, ina, matom 
  integer :: alpha_iatomha2, alpha2, nssh_tot2
  integer, dimension (:), allocatable :: mapindex
  integer, dimension (:), allocatable :: B
  integer, dimension (:,:), allocatable :: M
 
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
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

  alpha_iatomha2 = 0 
  
  print*,'alpha,alpha_iatomha2'
  do alpha = 1, nssh_tot
    if (fix_shell_charge(alpha) .eq. 0) alpha_iatomha2=alpha_iatomha2+1
    mapindex(alpha)=alpha_iatomha2
    print*,alpha,alpha_iatomha2
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
                      & + gshell_hh(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & + gshell_xc(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & + gshell_xc(alpha, beta  ,alpha_iatom, matom, beta_iatom) &
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
                      & + gshell_hh(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & + gshell_xc(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & + gshell_xc(alpha, beta  ,alpha_iatom, matom, beta_iatom) &
                      & - fshell_xc(beta , alpha ,alpha_iatom, matom, beta_iatom) &
                      & - fshell_xc(alpha, beta  ,alpha_iatom, matom, beta_iatom))&
                     & + exc_shell_aa(alpha)- vxc_shell_aa(alpha)  
          
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
    real function gshell_hh(alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom)
      use M_system, only: g_h, get_orb_ofshell, get_l_ofshell, get_issh_ofshell
      implicit none
      integer, intent(in) :: alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom
      integer :: mu_min, mu_max, imu_local
      gshell_hh = 0.0
      mu_min = get_orb_ofshell(alpha)
      mu_max = mu_min + 2*get_l_ofshell(alpha)
      do imu_local = mu_min, mu_max
        gshell_hh = gshell_hh + g_h(imu_local, imu_local, get_issh_ofshell(beta), alpha_iatom, neigh_beta_iatom, beta_iatom)
      end do
      gshell_hh = gshell_hh /  (2*get_l_ofshell(alpha) + 1)
    end function gshell_hh
 
    real function gshell_xc(alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom)
      use M_system, only: g_xc, get_orb_ofshell, get_l_ofshell, get_issh_ofshell
      implicit none
      integer, intent(in) :: alpha, beta, alpha_iatom, neigh_beta_iatom, beta_iatom
      integer :: mu_min, mu_max, imu_local
      gshell_xc = 0.0
      mu_min = get_orb_ofshell(alpha)
      mu_max = mu_min + 2*get_l_ofshell(alpha)
      do imu_local = mu_min, mu_max
        gshell_xc = gshell_xc + g_xc(imu_local, imu_local, get_issh_ofshell(beta), alpha_iatom, neigh_beta_iatom, beta_iatom)
      end do
      gshell_xc = gshell_xc /  (2*get_l_ofshell(alpha) + 1)
    end function gshell_xc

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




end subroutine stationary_charges
