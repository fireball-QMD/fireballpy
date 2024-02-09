subroutine neighbors ()
  use M_system
  use M_fdata, only: nssh, rcutoff
  implicit none
  integer :: iatom,jatom,mbeta,num_neigh,num_neigh_vdw,in1,imu,in2,neighcount
	real :: distance2,rcutoff_j, rcutoff_i,distance,range2,rc_max
  integer :: mbeta_max2
  integer :: neigh_max_old, ii, jj
  if (icluster .eq. 1) mbeta_max = 0
  rc_max = 0.00
  do iatom = 1, natoms
    num_neigh = 0
    rcutoff_i = 0.0d0
    in1 = imass(iatom)
    do imu = 1, nssh(in1)
      if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do
    do mbeta = 0, mbeta_max
      do jatom = 1, natoms
        rcutoff_j = 0.0d0
        in2 = imass(jatom)
        do imu = 1, nssh(in2)
          if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
        end do
        distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2+ (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
        distance = sqrt(distance2)
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
        rc_max = max(rc_max,(rcutoff_i + rcutoff_j))
        if (distance2 .le. range2) then
          if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and.iatom .ne. jatom) then
            write (*,*) ' WARNING - atoms dangerously close! '
            write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
          end if
          num_neigh = num_neigh + 1
          neigh_j(num_neigh,iatom) = jatom
          neigh_b(num_neigh,iatom) = mbeta
        end if
      end do
   end do
   neighn(iatom) = num_neigh
 end do
end subroutine

