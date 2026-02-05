subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qin, Qout, &
  & rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, get_l_ofshell, &
  & get_orb_ofshell, get_issh_ofshell, g_h, g_xc, get_iatom_ofshell, ztot, &
  & g_h_shell, g_xc_shell,f_xc_shell,exc_aa_shell,vxc_aa_shell
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer beta_iatom            
  integer imu, inu          
  integer in1, in2          
  integer issh, jssh
  integer ineigh ,jatom         
  integer mbeta
  real(double),dimension(nssh_tot,nssh_tot) :: A
  real(double),dimension(nssh_tot) :: c, SQ ! carga
  real(double) :: Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real(double) aux
  integer :: beta, alpha, alpha_iatom, ina, matom 
  integer :: alpha_iatom2, alpha2, nssh_tot2, lwork, info
  integer, dimension (:), allocatable :: mapindex, ipiv
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
          M(mapindex(alpha),mapindex(beta)) = M(mapindex(alpha),mapindex(beta)) & 
          & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
          & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  )!&
          do ineigh = 1, neighn(beta_iatom)
            mbeta = neigh_b(ineigh,beta_iatom)
            jatom = neigh_j(ineigh,beta_iatom)
            in2 = imass(jatom)
            in1 = imass(beta_iatom)
            do imu = 1, num_orb(in1)
              do inu = 1, num_orb(in2)
                aux = g_h(imu,inu,issh,alpha_iatom,ineigh,beta_iatom) + g_xc(imu,inu,issh,alpha_iatom,ineigh,beta_iatom)
                B(mapindex(alpha)) = B(mapindex(alpha)) + rho(imu,inu,ineigh,beta_iatom)*aux
              end do ! inu
            end do ! imu
          end do ! ineigh
        endif !fix_shell_charge(beta) = 0
        if (fix_shell_charge(beta) .eq. 1) then
          B(mapindex(alpha)) = B(mapindex(alpha))- Qin(get_issh_ofshell(beta),beta_iatom)*( &
          & +  g_h_shell(beta , alpha ) + g_xc_shell(beta , alpha ) + g_xc_shell(alpha, beta  ) &
          & - f_xc_shell(beta , alpha ) - f_xc_shell(alpha, beta  ))!&
        endif !fix_shell_charge(beta) = 1
      end do !beta
      B(mapindex(alpha)) = B(mapindex(alpha)) + exc_aa_shell(alpha) - vxc_aa_shell(alpha)
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
    allocate(work(1), ipiv(nssh_tot))
    call dsysv( 'U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, -1, info )
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsysv( 'U', nssh_tot, 1, M, nssh_tot, ipiv, B, nssh_tot, work, lwork, info )
    deallocate(work)
    call dgesv(nssh_tot,1,M,nssh_tot,ipiv,B,nssh_tot,info )
    deallocate(ipiv)
  contains

    subroutine load_M() !Mx=B
      use M_system, only: g_xc_shell,f_xc_shell,imass,vxc_aa_shell,exc_aa_shell,nssh_tot, &
      & g_h_shell, natoms, get_shell_ofatom_imu, get_shell_ofatom_issh
      use M_fdata, only: gxc_1c,fxc_1c,exc_1c_0,vxc_1c_0, nssh
      implicit none
      integer :: iatom,count ,issh,kssh,alpha,beta,imu, matom, katom 
      f_xc_shell=0.0d0
      count=0
      do iatom=1,natoms
        do issh=1,nssh(imass(iatom))
          alpha = count + issh          
          do kssh=1, nssh(imass(iatom))
            beta = count + kssh
            f_xc_shell(alpha,beta) = fxc_1c(imass(iatom),issh,kssh)
          end do
          exc_aa_shell(alpha) = exc_1c_0(imass(iatom),issh)
        end do
        count = count + nssh(imass(iatom))
      end do

      g_h_shell = 0.0
      g_xc_shell=0.0d0
      vxc_aa_shell = 0.0d0
      do iatom = 1, natoms
        do imu=1,num_orb(imass(iatom))
          alpha = get_shell_ofatom_imu(iatom,imu) 
          matom=neigh_self(iatom) 
          do katom = 1, natoms
            do kssh = 1, nssh(imass(katom))
              beta = get_shell_ofatom_issh(katom,kssh)
              g_h_shell(alpha,beta) = g_h_shell(alpha,beta) +&
              & g_h(imu,imu, kssh, katom, matom, iatom) / (2*get_l_ofshell(alpha) + 1)
              g_xc_shell(alpha,beta) = gxc_1c(imass(iatom),imu,issh,kssh) / (2*get_l_ofshell(alpha) + 1)
            end do
          end do
          vxc_aa_shell(alpha) = vxc_aa_shell(alpha) + vxc_1c_0(imass(iatom),imu,imu) / (2*get_l_ofshell(alpha) + 1)
        end do
      end do
    end subroutine  load_M

end subroutine stationary_charges
