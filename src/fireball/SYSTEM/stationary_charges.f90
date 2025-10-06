subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qout, rho, nssh_tot, neigh_self,neigh_b, fix_shell_charge, getlimu, getimu
  use M_fdata, only: num_orb,nssh,lssh
  implicit none
  integer iatom            
  integer imu, inu          
  integer in1, in2          
  integer issh, jssh
  integer ineigh ,jatom         
  integer mbeta
  real,dimension(nssh_tot,nssh_tot) :: A
  real,dimension(nssh_tot) :: c, SQ ! carga
  real :: diff_err,Ep2
  integer issh1, mu_min, mu_max, l, inumorb
  real auxgS
  integer :: beta, alpha, ialp, ina, matom
  SQ = 0.0d0
  c = 0.0d0
  alpha = 0
!Borrar es solo para hacer test----------
  !do issh=1,nsh_max
  !  fix_shell_charge(issh)=fix_shell_charge_aux(issh)
  !end do
  !lo usamos para H2O HsHsOsp
  !                   s s s p 
  fix_shell_charge = [0,0,1,0]
!-----------------------------------------
! Mx=B x carga en shells q^mu
! B^mu = sum_ij [ P^mu_ij (g_h^mu_ij +g_xc^mu_ij) + exc_aa^mu -  vxc_aamu
! M_mu_i = D^mu_ii = g_h^i_mumu + g_xc^mu_ii + g_xc^i_mumu - f_xc^mu_ii - f_xc^i_mumu
  do ialp = 1, natoms
    ina = imass(ialp)
    do issh = 1, nssh(ina)
      alpha = alpha + 1 ! map orbitales -> shells
      beta = 0
      do iatom = 1, natoms
        in1 = imass(iatom)
        matom = neigh_self(iatom)
        do issh1 = 1, nssh(in1)
          beta = beta + 1
          A(alpha,beta) = A(alpha,beta) & 
                      & + gshell_hh(beta, issh,ialp,matom,iatom) 


! g_h(imu,imu,issh,ialp,matom,iatom) &
!            &               + g_xc(imu,imu,issh,ialp,matom,iatom) &
!            &               + g_

            !g_h^i_mumu + g_xc^mu_ii + g_xc^i_mumu - f_xc^mu_ii - f_xc^i_mumu
            ! antes sumamba :  gvhxc(imu,imu,issh,ialp,matom,iatom)
            !M(alpha,beta) =  auxgS
         ! end do ! end do issh1
          do ineigh = 1, neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            do imu = 1, num_orb(in1)
              do inu = 1, num_orb(in2)
                c(alpha) = c(alpha) + 0.0d0 !&
!             &      rho(imu,inu,ineigh,iatom)*gvhxc(imu,inu,issh,ialp,ineigh,iatom)
              end do ! end do inu
            end do ! end do imu
          end do ! end do ineigh
        end do ! end do issh1
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
  stop    
contains
  real function gshell_hh(alpha, issh, ialp, matom, iatom)
    use M_system, only: g_h, getimu, getlimu
    implicit none
    integer, intent(in) :: alpha, issh, ialp, matom, iatom
    integer :: mu_min, mu_max, imu_local
    gshell_hh = 0.0
    mu_min = getimu(alpha)
    mu_max = mu_min + 2*getlimu(alpha)
  
  ! Sumar sobre los orbitales
  do imu_local = mu_min, mu_max
    gshell_hh = gshell_hh + g_h(imu_local, imu_local, issh, ialp, matom, iatom)
  end do
  gshell_hh = gshell_hh /  (2*getlimu(alpha) + 1)
  
end function gshell_hh
 
end subroutine stationary_charges
