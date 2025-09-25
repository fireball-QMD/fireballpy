subroutine stationary_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_j, neighn, numorb_max, Qout, rho, nssh_tot, neigh_self,neigh_b
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
  do ialp = 1, natoms
    ina = imass(ialp)
    do issh = 1, nssh(ina)
      alpha = alpha + 1 ! map orbitales -> shells
      beta = 0
      do iatom = 1, natoms
        in1 = imass(iatom)
        matom = neigh_self(iatom)
        inumorb = 1
        do issh1 = 1, nssh(in1)
          beta = beta + 1
          l = lssh(issh1,in1)
          mu_min = inumorb
          mu_max = mu_min+2*l
          do imu = mu_min, mu_max
            auxgS =  auxgS  + 0.0 !  gvhxc(imu,imu,issh,ialp,matom,iatom)
            auxgS = auxgS/(2*l+1) 
            !M(alpha,beta) =  auxgS
            A(alpha,beta) = auxgS
            inumorb = inumorb + 2*l+1
          end do ! end do issh1
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
end subroutine stationary_charges
