!    This routine takes the input density and determines the
! exchange-correlation potential and energy atom-like matrix elements 
! in the average-density approx. (see PRB 40, 3979 (1989) ).
! The double-counting xc-contribution is also calculated here
! Formula:
! <i,mu|V_xc(n)|i,nu> = V_xc(n^a)<i,mu|i,nu> + V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>] -
!   - <i,mu|V_xc(n_i)|i,nu> + V_xc(n_i^a)<i,mu|i,nu> +
!   + V'_xc(n_i^2)[<i,mu|n_i|i,nu> + n_i^a<i,mu|i,nu>] 
! where
!  V'_xc = d(V_xc)/dn
!  n^a ... average density
!  n_i ... density on i-site
! 
subroutine build_ca_olsxc_on (in1, iatom, bcxcx, xc)
  use M_system
  implicit none
  integer, intent (in) :: in1
  integer, intent (in) :: iatom
  real, intent (out), dimension (numorb_max, numorb_max) :: bcxcx
  real, intent (out) :: xc

  integer imu
  integer ind1
  integer ind2
  integer inu
  integer issh
  integer jssh
  integer l1
  integer l2
  integer n1
  integer n2

  real dexc
  real d2exc 
  real dmuxc
  real d2muxc
  real exc
  real muxc
  real dexci
  real d2exci
  real dmuxci
  real d2muxci
  real exci
  real muxci
  real q_mu

  real, dimension (nsh_max,nsh_max) :: arho
  real, dimension (nsh_max,nsh_max) :: arhoi
  real, dimension (numorb_max, numorb_max) :: denx
  real, dimension (numorb_max, numorb_max) :: deni

  xc = 0.0d0
  bcxcx = 0.0d0
  do inu = 1, nssh(in1)
    do imu = 1, nssh(in1)
      arho(imu,inu) = arho_on(imu,inu,iatom)
      arhoi(imu,inu) = arhoi_on(imu,inu,iatom)
    end do
  end do
  do inu = 1, num_orb(in1)
    do imu = 1, num_orb(in1)
      denx(imu,inu) = rho_on(imu,inu,iatom) 
      deni(imu,inu) = rhoi_on(imu,inu,iatom)
    end do
  end do
  n1 = 0
  do issh = 1, nssh(in1)
    l1 = lssh(issh,in1)
    n1 = n1 + l1 + 1
    call cepal (arhoi(issh,issh), exci, muxci, dexci, d2exci, dmuxci,  d2muxci)
    call cepal (arho(issh,issh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
    ! Now calculate different parts of the XC-matrix 
    ! Formula:
    !  <i,mu|V_xc(n)|i,nu> =
    !  V_xc(n^a)<i,mu|i,nu> + V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>] -
    !   - <i,mu|V_xc(n_i)|i,nu> + V_xc(n_i^a)<i,mu|i,nu> +
    !   + V'_xc(n_i^2)[<i,mu|n_i|i,nu> - n_i^a<i,mu|i,nu>] 
    !  where
    !  V_xc(n^a)    ... muxc
    !  V'_xc(n^a)   ... dmuxc
    !  V_xc(n_i^a)  ... muxc0
    !  V'_xc(n_i^a) ... dmuxc0
    !
    !  n^a    ... arho
    !  n_i    ... arhoi
    !
    !  mu /= nu :   <mu|nu> = 0 (off-diagonal term) 
    !  mu = nu  :   <mu|nu> = 1 (diagonal term)
    !
    !  keeping in mind that <i,mu|V_xc(n_i)|i,nu> is done in assemble_1c()   !  !  !  
    do ind1 = -l1, l1
      imu = n1 + ind1
      ! V_xc(n^a)<i,mu|i,nu>
      bcxcx(imu,imu) = muxc
      ! V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>] 
      bcxcx(imu,imu) =  bcxcx(imu,imu) + dmuxc*(denx(imu,imu) - arho(issh,issh) )
      ! - V_xc(n_i^a)<i,mu|i,nu> -  V'_xc(n_i^2)[<i,mu|n_i|i,nu> - n_i^a<i,mu|i,nu>]}
      bcxcx(imu,imu) = bcxcx(imu,imu) - muxci - dmuxci*(deni(imu,imu) - arhoi(issh,issh) )
    end do
    ! Contribution to the double-counting correction
    ! average charge per shell
    q_mu = Qin(issh,iatom) / (2*l1 + 1)
    do ind1 = -l1, l1
      imu = n1 + ind1
      ! SNXC part
      xc = xc + q_mu*(exc - muxc + (dexc - dmuxc)*(denx(imu,imu) - arho(issh,issh)))
      ! OLSXC part
      xc = xc + q_mu*(muxci - exci - (dexci-dmuxci)*(deni(imu,imu)-arhoi(issh,issh)) )
    end do
    n1 = n1 + l1
  end do   ! end 'do issh = 1, nssh(in1)'
  ! END - DIAGONAL TERMS

  ! OFF-DIAGONAL TERMS
  ! Calculate the vxc-matrix elements for the "non-diagonal" terms
  n1 = 0
  do issh = 1, nssh(in1)
    l1 = lssh(issh,in1)
    n1 = n1 + l1 + 1
    n2 = 0
    do jssh = 1, nssh(in1)
      l2 = lssh(jssh,in1)
      n2 = n2 + l2 + 1
      call cepal(arhoi(issh,jssh), exci, muxci, dexci, d2exci, dmuxci, d2muxci)
      call cepal(arho(issh,jssh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
      do ind1 = -l1, l1
        imu = n1 + ind1
        do ind2 = -l2, l2
          inu = n2 + ind2
          if (imu .ne. inu) then
            ! V_xc(n^a)<i,mu|i,nu> = 0
            ! V'_xc(n^a) n^a <i,mu|i,nu> = 0
            ! V_xc(n_i^a)<i,mu|i,nu> = 0
            ! V'_xc(n_i^2) n_i^a <i,mu|i,nu> = 0
            ! V'_xc(n^a) <i,mu|n|i,nu>
            bcxcx(imu,inu) = dmuxc*denx(imu,inu) 
            ! - V'_xc(n_i^2)<i,mu|n_i|i,nu> 
            bcxcx(imu,inu) = bcxcx(imu,inu) - dmuxci*deni(imu,inu) 
          end if
        end do   !do ind2 = -l2, l2
      end do   !do ind1 = -l1, l1
      n2 = n2 + l2
    end do   !do jssh = 1, nssh(in1)
    n1 = n1 + l1
  end do   !do issh = 1, nssh(in1)
  return
end subroutine build_ca_olsxc_on
