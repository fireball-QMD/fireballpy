subroutine neighbors ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: icluster, natoms, ratom, imass, mbeta_max, neigh_b, neigh_j, neighn, xl
  use M_fdata, only: nssh, rcutoff
  implicit none
  integer :: iatom,jatom,mbeta,num_neigh,in1,imu,in2
  real(double) :: distance2,rcutoff_j, rcutoff_i,distance,range2,rc_max

  if (icluster .eq. 1) mbeta_max = 0
  rc_max = 0.00d0
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
        distance2 = sum((ratom(:,jatom) + xl(:,mbeta) - ratom(:,iatom))**2)
        distance = sqrt(distance2)
        range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
        rc_max = max(rc_max,(rcutoff_i + rcutoff_j))
        if (distance2 .le. range2) then
          num_neigh = num_neigh + 1
          neigh_j(num_neigh,iatom) = jatom
          neigh_b(num_neigh,iatom) = mbeta
        end if
      end do
   end do
   neighn(iatom) = num_neigh
 end do
end subroutine neighbors
