subroutine common_neighborsPP ()
  use M_system
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer imu
  integer in1, in2, indna
  integer ineigh
  integer jatom
  integer jbeta
  integer jneigh
  integer katom
  integer kbeta
  integer kneigh
  integer num_neigh
  integer mbeta

  real*8  distance
  real*8, dimension (3) :: dvec
  real*8, dimension (3) :: diff
  real*8, dimension (3) :: vec1, vec2, vec3, vec

  if(icluster .eq. 1) mbeta_max = 0
  do ialp = 1, natoms
    vec2(:) = ratom(:,ialp)
    num_neigh = 0
    do ineigh = 1, nPPxn(ialp)
      iatom = nPPx_j(ineigh,ialp)
      ibeta = nPPx_b(ineigh,ialp)
      in1 = imass(iatom)
      if (.not. (iatom .eq. ialp .and. ibeta .eq. 0)) then
        vec1(:) = ratom(:,iatom) + xl(:,ibeta)
        do jneigh = 1, nPPxn(ialp)
          jatom = nPPx_j(jneigh,ialp)
          jbeta = nPPx_b(jneigh,ialp)
          in2 = imass(jatom)
          if (.not. (jatom .eq. ialp .and. jbeta .eq. 0) .and. (ineigh .ne. jneigh)) then
            vec3(:) = ratom(:,jatom) + xl(:,jbeta)
            num_neigh = num_neigh + 1
            if (num_neigh .gt. neighPP_max**2) then
              write (*,*) ' num_neigh = ', num_neigh, ' within neighbors.f90! '
              stop
            end if ! if(num_neigh)
            neighPP_comj(1,num_neigh,ialp) = iatom
            neighPP_comb(1,num_neigh,ialp) = ibeta
            neighPP_comj(2,num_neigh,ialp) = jatom
            neighPP_comb(2,num_neigh,ialp) = jbeta
            ! We also need to know for a given ialp and (iatom,ibeta), what is the m value
            ! for (jatom,jbeta) with respect to iatom. That is, jatom is the m'th neighbor
            ! of iatom. What is m?
            diff = vec3 - vec1
            neighPP_comm(num_neigh,ialp) = -9999
            do kneigh = 1, neighPPn(iatom)
              katom = neighPP_j(kneigh,iatom)
              kbeta = neighPP_b(kneigh,iatom)
              vec(:) = xl(:,kbeta) + ratom(:,katom) - ratom(:,iatom)
              if ((abs(vec(1) - diff(1)) .lt. 1.0d-4) .and. (abs(vec(2) - diff(2)) .lt. 1.0d-4) .and.(abs(vec(3) - diff(3)) .lt. 1.0d-4))  neighPP_comm(num_neigh,ialp) = kneigh
            end do ! do kneigh
            if (neighPP_comm(num_neigh,ialp) .eq. -9999) num_neigh = num_neigh - 1 
          end if ! if(.not.(jatom .eq. ialp)
        end do ! do jneigh 
      end if ! if(.not.(iatom .eq. ialp) 
    end do ! do ineigh 
    neighPP_comn(ialp) = num_neigh
  end do


  return
end subroutine common_neighborsPP
