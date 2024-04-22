subroutine MULLIKEN_CHARGES()        
 use M_system
 use M_constants
 use M_fdata, only: num_orb,nssh,lssh
 implicit none
 integer iatom            
 integer ikpoint          
 integer imu, inu          
 integer in1, in2          
 integer issh, jssh
 integer ineigh ,jatom,jneigh          
 integer noccupy  
 integer mqn             
 real(8), dimension (numorb_max, natoms) :: QMulliken
 QMulliken = 0.0d0             
 do iatom = 1, natoms
   in1 = imass(iatom)
   do ineigh = 1, neighn(iatom)
     jatom = neigh_j(ineigh,iatom)
     in2 = imass(jatom)
     jneigh = neigh_back(iatom,ineigh)
     do imu = 1, num_orb(in1)
      do inu = 1, num_orb(in2)
       QMulliken(imu,iatom) = QMulliken(imu,iatom)+ 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
      end do
     end do
   end do 
   imu = 0
   do issh = 1, nssh(in1)
     do mqn = 1, 2*lssh(issh,in1) + 1
      imu = imu + 1
      Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
     end do
     QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(issh,iatom)
   end do
 end do  
end subroutine MULLIKEN_CHARGES 
