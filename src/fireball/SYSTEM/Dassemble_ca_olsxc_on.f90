subroutine Dassemble_ca_olsxc_on ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, imass, neigh_b, neigh_j, neighn, neigh_self, numorb_max, Qin, rho, rho_on, arho_on, arhop_on, rhop_on, &
    & faxc, faxc_ca, dxcdcc
  use M_fdata, only: nssh, lssh
  implicit none
  integer iatom
  integer imu
  integer in1
  integer in2
  integer index
  integer ineigh
  integer inu
  integer issh
  integer ix
  integer jatom
  integer jndex
  integer jssh
  integer l1 
  integer l2
  integer matom
  integer mbeta
  integer n1
  integer n2
  real(double) exc
  real(double) dexc
  real(double) d2exc
  real(double) d2muxc
  real(double) dmuxc
  real(double) muxc
  real(double) rhoxc
  real(double) rhoxc_av
  real(double) q_mu
  real(double), dimension (3, numorb_max, numorb_max) :: bcxcpx
  real(double), dimension (numorb_max, numorb_max) :: delta
  real(double), dimension (3) :: drhoxc
  real(double), dimension (3) :: drhoxc_av
  faxc = 0.0d0
  faxc_ca = 0.0d0
  dxcdcc  = 0.0d0
  delta = 0.0d0
  do imu = 1,numorb_max
    delta(imu,imu) = 1.0d0
  enddo
  do iatom = 1, natoms
    in1 = imass(iatom)
    matom = neigh_self(iatom)
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      in2 = imass(jatom)
      bcxcpx = 0.0d0
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
      else
        n1 = 0
        do issh = 1, nssh(in1)
          l1 = lssh(issh,in1)
          n1 = n1 + l1 + 1
          rhoxc_av = arho_on(issh,issh,iatom)
          drhoxc_av(:) = arhop_on(:,issh,issh,ineigh,iatom)
          call cepal(rhoxc_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
          q_mu = Qin(issh,iatom) / (2.0d0*l1 + 1)
          do index = -l1, l1
            imu = n1 + index
            rhoxc  = rho_on(imu,imu,iatom)
            drhoxc(:) = rhop_on(:,imu,imu,ineigh,iatom)
            bcxcpx(:,imu,imu) =   ((dexc - dmuxc)*drhoxc_av(:) + (d2exc - d2muxc)*(rhoxc - rhoxc_av)*drhoxc_av(:)  +  (dexc - dmuxc)*(drhoxc(:) - drhoxc_av(:)) )*q_mu
            dxcdcc(:,ineigh,iatom) =  dxcdcc(:,ineigh,iatom) -  bcxcpx(:,imu,imu)
          end do !do index = -l1, l1
          n1 = n1 + l1
        end do ! end do issh
        n1 = 0
        do issh = 1, nssh(in1)
          l1 = lssh(issh,in1)
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, nssh(in1)
            l2 = lssh(jssh,in1)
            n2 = n2 + l2 + 1
            ! calculate XC potentials 
            rhoxc_av      = arho_on(issh,jssh,iatom)
            drhoxc_av(:) = arhop_on(:,issh,jssh,ineigh,iatom) 
            call cepal(rhoxc_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
            do index = -l1, l1
              imu = n1 + index
              do jndex = -l2, l2
                inu = n2 + jndex
                rhoxc  = rho_on(imu,inu,iatom)
                drhoxc(:) = rhop_on(:,imu,inu,ineigh,iatom)
                do ix = 1,3
                  bcxcpx(ix,imu,inu) = drhoxc(ix)*dmuxc + rhoxc*drhoxc_av(ix)*d2muxc  - delta(imu,inu)*drhoxc_av(ix)*d2muxc*rhoxc_av                 
                  faxc_ca(ix,ineigh,iatom) = faxc_ca(ix,ineigh,iatom)  -  bcxcpx(ix,imu,inu)*rho(imu,inu,matom,iatom)
                end do 
              end do 
            end do 
            n2 = n2 + l2
          end do 
          n1 = n1 + l1
        end do 
      end if
    end do !  ineigh
  end do !  iatom

end subroutine Dassemble_ca_olsxc_on
