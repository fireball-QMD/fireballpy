subroutine Dassemble_2c_PP ()
  use M_system
  use M_fdata, only: num_orb, num_orbPP
  implicit none
  integer iatom
  integer ierror
  integer imu
  integer in1
  integer in2
  integer ineigh
  integer inu
  integer ix
  integer jatom
  integer jneigh
  integer kneigh
  integer matom
  integer mbeta
  integer mneigh_self
  integer ncc
  real(8), dimension (numorb_max) :: cl
  real(8), dimension (3, numorb_max, numorb_max) :: fnlb
  fanl = 0.0d0
  fotnl = 0.0d0
  do iatom = 1, natoms
    matom = neighPP_self(iatom)
    in1 = imass(iatom)
    do ineigh = 1, nPPn(iatom)    
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
      in2 = imass(jatom)
      call cl_value(in2,cl)
      mneigh_self = nPP_self(iatom)
      do inu = 1, num_orb(in1)
        do imu = 1, num_orb(in1)
          fnlb(:,imu,inu) = 0.0d0
          if (ineigh .ne. mneigh_self) then
            do ncc = 1, num_orbPP(in2)
              do ix = 1,3
                fnlb(ix,imu,inu) = fnlb(ix,imu,inu)  + cl(ncc)*(spVNL(ix,imu,ncc,ineigh,iatom)*sVNL(inu,ncc,ineigh,iatom)  + sVNL(imu,ncc,ineigh,iatom)*spVNL(ix,inu,ncc,ineigh,iatom))
              end do 
            end do 
          end if 
        end do 
      end do 
      do inu = 1, num_orb(in1)
        do imu = 1, num_orb(in1)
          do ix = 1,3
            fanl(ix,ineigh,iatom) =  fanl(ix,ineigh,iatom) - rhoPP(imu,inu,matom,iatom)*fnlb(ix,imu,inu)
          end do
        end do 
      end do 
    enddo 

    do ineigh = 1, nPPxn(iatom)
      mbeta = nPPx_b(ineigh,iatom)
      jatom = nPPx_j(ineigh,iatom)
      in2 = imass(jatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
      else
        call cl_value(in1,cl)
        mneigh_self = nPP_self(iatom)
        jneigh = nPPx_point(ineigh,iatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            fnlb(:,imu,inu) = 0.0d0
            do ncc = 1, num_orbPP(in1)
              do ix = 1,3
                fnlb(ix,imu,inu) = fnlb(ix,imu,inu)  - cl(ncc)*sVNL(imu,ncc,mneigh_self,iatom)*spVNL(ix,inu,ncc,jneigh,jatom)
              end do 
            end do 
          end do 
        end do 
        kneigh = nPPx_map(ineigh,iatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            do ix = 1,3
              fotnl(ix,ineigh,iatom) = fotnl(ix,ineigh,iatom) - rhoPP(imu,inu,kneigh,iatom)*fnlb(ix,imu,inu)
            end do 
          end do 
        end do 
      end if 
    end do 
  end do  
end subroutine Dassemble_2c_PP

