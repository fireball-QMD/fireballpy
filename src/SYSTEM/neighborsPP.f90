!  Find all neighbors to atoms in the central cell in case of Kleiman Bylander 
! fully separable pseudopotential. The distance of neighbor PP pair must be 
! less then sum of wave function and PP radius cutoff   
! An atom is at lattice vector xl(mbeta), and basis ratom(iatom).  
! We refer to this as (mbeta,iatom). An atom in the central cell is 
! at (0,iatom).
! Find all PP-neighbors to atom (0,iatom).
!
! neighnPP(i)=# of neighbors of atom i
! neighjPP(i,m)= j-sub-m, the j value of the m'th neighbor.
! neighbPP(i,m)= beta-sub-m, the beta value for the m'th neighbor.
!
!       The important quantity here is neighj, which indicates which basis
! vector we have: This identifies the species a through the array 
! imass (neighjPP(iatom,ineigh)) is the type (1 or 2) of the ineigh'th 
! PP-neighbor u of atom iatom.  Furthermore, this atom is located at 
! xl(:,neighbPP(iatom,ineigh)) + ratom(:,neighjPP(iatom,ineigh)).

subroutine neighborsPP ()
  use M_system
  use M_fdata,only : nssh, rcutoff, rc_PP    
  implicit none
  integer iatom
  integer ineigh
  integer ibeta
  integer imu
  integer in1
  integer in2
  integer jatom
  integer jneigh
  integer jbeta
  integer katom
  integer kneigh
  integer kbeta
  integer ialp
  integer mbeta
  integer num_neigh
  real(8) distance
  real(8) distance2
  real(8) range2
  real(8) rcutoff_i
  real(8) rcutoff_j
  real(8), dimension (3) :: vec1, vec2, vec3, vec
  logical flag

  neighPP_j = 0.0d0
  neighPP_b  = 0.0d0

  !   LIST nPP
  ! First we'are going to build nPP list of neighbors. The list includes all
  ! neighbors of iatom with nonzero overlap <phi_i|Psi_j>, where
  ! phi_i .. is atomic wave function and Psi_j is pseudo wave function
  if (icluster .eq. 1) mbeta_max = 0
  do iatom = 1, natoms
   num_neigh = 0
   rcutoff_i = 0.0d0
   in1 = imass(iatom)
   do imu = 1, nssh(in1)
    if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
   end do
   do mbeta = 0, mbeta_max
    do jatom = 1, natoms
     in2 = imass(jatom) 
     rcutoff_j = rc_PP(in2)
     distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &    + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &    + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
     distance = sqrt(distance2)
     range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
     if (distance2 .le. range2) then
      if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and. iatom .ne. jatom) then
       write (*,*) ' WARNING - atoms dangerously close! '
       write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
      end if 
      num_neigh = num_neigh + 1
      nPP_j(num_neigh,iatom) = jatom
      nPP_b(num_neigh,iatom) = mbeta
     end if ! if(distance2 .le. range2)
    end do ! do jatom 
   end do ! do mbeta
   nPPn(iatom) = num_neigh
  end do ! do iatom

  ! LIST nPPx
  ! These two lists are not identical if we have more species in the system vwith different Rc and Rc_PP
  do iatom = 1, natoms
   num_neigh = 0
   in1 = imass(iatom)
   rcutoff_i = rc_PP(in1)
   do mbeta = 0, mbeta_max
    do jatom = 1, natoms
     rcutoff_j = 0.0d0
     in2 = imass(jatom) 
     do imu = 1, nssh(in2)
      if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
     end do
     distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &    + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &    + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
     distance = sqrt(distance2)
     range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
     if (distance2 .le. range2) then
      if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and. iatom .ne. jatom) then
       write (*,*) ' WARNING - atoms dangerously close! '
       write (*,*) ' WARNING - atoms dangerously close! '
       write (*,*) ' WARNING - atoms dangerously close! '
       write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
      end if 
      num_neigh = num_neigh + 1
      nPPx_j(num_neigh,iatom) = jatom
      nPPx_b(num_neigh,iatom) = mbeta
     end if ! if(distance2 .le. range2)
    end do ! do jatom 
   end do ! do mbeta
   nPPxn(iatom) = num_neigh
  end do ! do iatom

 !   LIST neighPP
  do iatom = 1, natoms
   num_neigh = 0
   vec1(:) = ratom(:,iatom)
   do ineigh = 1, nPPn(iatom)
     ialp = nPP_j(ineigh,iatom)
     ibeta = nPP_b(ineigh,iatom)
     vec2(:) = ratom(:,ialp) + xl(:,ibeta)
     do jneigh = 1, nPPxn(ialp)
       jatom = nPPx_j(jneigh,ialp)
       jbeta = nPPx_b(jneigh,ialp)
       vec3(:) = xl(:,jbeta) + xl(:,ibeta)
       do mbeta = 0 , mbeta_max
         vec(:) = xl(:,mbeta)
         distance = sqrt( ( vec(1) - vec3(1) )**2.0d0  +   &
                &         ( vec(2) - vec3(2) )**2.0d0  +   &
                &         ( vec(3) - vec3(3) )**2.0d0  ) 
         if(distance .lt. 0.0001) exit
       enddo
       if (mbeta .gt. mbeta_max) then 
         write (*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++'
         write (*,*) ' mbeta_max =', mbeta_max
         write (*,*) ' vec_look =', vec3(:)
         write (*,*) ' Fireball cannot find desired neighboring cell. '
         write (*,*) ' This propably means you need to increase number '
         write (*,*) ' of the periodic unit cell box and to recompile the code.'
         write (*,*) ' For details see INITIALIZERS/initboxes.f90 file.'
         write (*,*) ' Fireball is going to die!'
         write (*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++'
         stop
       endif
       flag = .false.  
       do kneigh = 1, num_neigh
         katom = neighPP_j(kneigh,iatom)
         kbeta = neighPP_b(kneigh,iatom)
         if( katom .eq. jatom .and. kbeta .eq. mbeta ) flag = .true.
       enddo ! do kneigh
       if(.not. flag ) then 
         num_neigh = num_neigh + 1
         if (num_neigh .gt. neighPP_max**2) then
           write (*,*) 'Oh no. In common_neighbors.f90, we have too many'
           write (*,*) 'common neighbors. In MODULES/dimensions.f90 '
           write (*,*) 'dimension neigh_max**2 = ', neighPP_max**2
           write (*,*) 'So far (*but still counting) we have '
           write (*,*) 'num_neigh = ',num_neigh,' within neighbors.f90!'
           stop
         end if
         neighPP_j(num_neigh,iatom) = jatom
         neighPP_b(num_neigh,iatom) = mbeta
       endif ! if(flag)
     end do ! do jneigh
   end do ! do ineigh
   neighPPn(iatom) = num_neigh
  end do ! do iatom
  do iatom = 1, natoms
    vec1(:) = ratom(:,iatom)
    do ineigh = 1,nPPxn(iatom)
       jatom = nPPx_j(ineigh,iatom)
       jbeta = nPPx_b(ineigh,iatom)
       do kneigh = 1,nPPn(jatom)
         katom = nPP_j(kneigh,jatom)
         kbeta = nPP_b(kneigh,jatom)
         vec2(:) = ratom(:,katom) + xl(:,kbeta) + xl(:,jbeta)
         distance = sqrt( ( vec1(1) - vec2(1) )**2.0d0  +      &
               &          ( vec1(2) - vec2(2) )**2.0d0  +      &
               &          ( vec1(3) - vec2(3) )**2.0d0  )
         if (distance .lt. 0.0001d0) then
           nPPx_point(ineigh,iatom) = kneigh
           exit
         endif ! if (distance)
       enddo ! do jneigh
    enddo ! do ineigh
  enddo ! do iatom

  !   MAP nPPx_map  (nPPx -> neighPP)
  do iatom = 1, natoms
     do ineigh = 1,nPPxn(iatom)
        jatom = nPPx_j(ineigh,iatom)
        jbeta = nPPx_b(ineigh,iatom)
        do kneigh = 1,neighPPn(iatom)
          katom = neighPP_j(kneigh,iatom)
          kbeta = neighPP_b(kneigh,iatom)
          if (katom .eq. jatom .and. kbeta .eq. jbeta) then
            nPPx_map(ineigh,iatom) = kneigh 
            exit
          endif ! if (katom)
        enddo ! do kneigh
     enddo ! do ineigh
  enddo ! do iatom
  !   MAP nPP_map  (nPP -> neighPP)
  do iatom = 1, natoms
     do ineigh = 1,nPPn(iatom)
        jatom = nPP_j(ineigh,iatom)
        jbeta = nPP_b(ineigh,iatom)
        do kneigh = 1,neighPPn(iatom)
          katom = neighPP_j(kneigh,iatom)
          kbeta = neighPP_b(kneigh,iatom)
          if(katom .eq. jatom .and. kbeta .eq. jbeta) then
            nPP_map(ineigh,iatom) = kneigh 
            exit
          endif ! if (katom)
        enddo ! do kneigh
     enddo ! do ineigh
  enddo ! do iatom

  !   SELF nPP_self 
  nPP_self = -999
  do iatom = 1, natoms
    do ineigh = 1, nPPn(iatom)
      mbeta = nPP_b(ineigh,iatom)
      jatom = nPP_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) nPP_self(iatom) = ineigh
    end do
  end do

  !   SELF nPP_self 
  nPPx_self = -999
  do iatom = 1, natoms
    do ineigh = 1, nPPxn(iatom)
      mbeta = nPPx_b(ineigh,iatom)
      jatom = nPPx_j(ineigh,iatom)
      if (iatom .eq. jatom .and. mbeta .eq. 0) nPPx_self(iatom) = ineigh
    end do
  end do
  return
end subroutine neighborsPP

