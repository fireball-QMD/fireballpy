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
  real(8) y
  real(8), dimension (numorb_max, natoms) :: QMulliken
  real(8), dimension (3) :: vec, r1, r2, r21
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
