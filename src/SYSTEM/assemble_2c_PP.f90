! This routine assembles all of the two-center and degenerate two-center interactions.
subroutine assemble_2c_PP ()
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
  integer isorp
  integer jatom
  integer jneigh
  integer kneigh
  integer matom
  integer mbeta
  integer mneigh_self
  integer ncc
 
  real*8, dimension (numorb_max) :: cl
  real*8, dimension (numorb_max, numorb_max) :: PPx

  vnl = 0.0d0

  do iatom = 1, natoms
    matom = neighPP_self(iatom)
    in1 = imass(iatom)
    do ineigh = 1, nPPn(iatom)    !  <==== loop over i's neighbors
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
      in2 = imass(jatom)
      call cl_value (in2, cl)
      ! in1 twice because it is an atom case.
      do inu = 1, num_orb(in1)
        do imu = 1, num_orb(in1)
          PPx(imu,inu) = 0.0d0
          do ncc = 1, num_orbPP(in2)
            PPx(imu,inu) = PPx(imu,inu)+ cl(ncc)*sVNL(imu,ncc,ineigh,iatom)*sVNL(inu,ncc,ineigh,iatom)
          end do
        end do
      end do
 
      ! Final assembly of vnl - the energy piece.
      do inu = 1, num_orb(in1)
        do imu = 1, num_orb(in1)
          vnl(imu,inu,matom,iatom) = vnl(imu,inu,matom,iatom) + PPx(imu,inu)
        end do
      end do
    enddo   ! do ineigh
  enddo   ! do iatom
  ! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
  do iatom = 1, natoms
    matom = neighPP_self(iatom)
    in1 = imass(iatom)
    do ineigh = 1, nPPxn(iatom) 
      mbeta = nPPx_b(ineigh,iatom)
      jatom = nPPx_j(ineigh,iatom)
      in2 = imass(jatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        if (nPPx_self(iatom) .ne. ineigh) then
          write (*,*) ' Something really wrong in assemble_2c_PP.f90 '
          write (*,*) ' iatom, jatom, mbeta = ', iatom, jatom, mbeta
          write (*,*) ' neigh_self(iatom), ineigh = ', nPPx_self(iatom), ineigh  
        end if   ! if(neighPP_self)
      else   ! if(iatom .eq. jatom)

        ! Case 1. PP is iatom.  <i | VNL(i) |j>.
        call cl_value (in1, cl)
        jneigh = nPPx_point(ineigh,iatom)
        ! <phi_i|Psi_i>  ->  nPP(mneigh_self,iatom)
        mneigh_self = nPP_self(iatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            PPx(imu,inu) = 0.0d0
            do ncc = 1, num_orbPP(in1)
              PPx(imu,inu) = PPx(imu,inu)+ cl(ncc)*sVNL(imu,ncc,mneigh_self,iatom)*sVNL(inu,ncc,jneigh,jatom)
            end do   ! do ncc
          end do   ! do imu
        end do   ! do inu
        ! Mapping to the global matrix
        kneigh = nPPx_map(ineigh,iatom)
        ! Assemble the global matrix 
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vnl(imu,inu,kneigh,iatom) = vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)
          end do   ! do imu
        end do   ! do inu
      endif   ! if(iatom .eq. jatom)
    enddo   ! do ineigh
  enddo   ! do iatom

  ! ASSEMBLE VNL ONTOP RIGHT CASE   <phi_i|Psi_j><Psi_j|phi_j>
  do iatom = 1, natoms
    matom = neighPP_self(iatom)
    in1 = imass(iatom)
    do ineigh = 1, nPPn(iatom) 
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
      in2 = imass(jatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
        if (nPP_self(iatom) .ne. ineigh) then
          write (*,*) ' Something real*8ly wrong in assemble_2c_PP.f90 '
          write (*,*) ' iatom, jatom, mbeta = ', iatom, jatom, mbeta
          write (*,*) ' neigh_self(iatom), ineigh = ',nPP_self(iatom), ineigh  
          stop
        end if   ! if(neighPP_self)
      else   ! if(iatom .eq. jatom)
        ! Now the second case. <i | V(j) | j>.
        call cl_value (in2, cl)
        ! Looking for <phi_j|Psi_j>, what is jneigh of jatom itself in the nPPx list 
        mneigh_self = nPP_self(jatom)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            PPx(imu,inu) = 0.0d0
            do ncc = 1, num_orbPP(in2)
              PPx(imu,inu) = PPx(imu,inu) + cl(ncc)*sVNL(imu,ncc,ineigh,iatom)*sVNL(inu,ncc,mneigh_self,jatom)
            end do
          end do
        end do

        kneigh = nPP_map(ineigh,iatom) 
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vnl(imu,inu,kneigh,iatom) = vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)
          end do
        end do
      end if   ! if(iatom .eq. jatom)
    end do   ! do ineigh
  end do   ! do iatom
  return
end subroutine assemble_2c_PP
