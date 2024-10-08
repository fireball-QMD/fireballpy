subroutine common_neighbors ()
  use iso_c_binding
  use M_system, only: natoms, ratom, imass, neigh_b, neigh_j, neighn, neigh_comb, neigh_comj, neigh_comm, neigh_comn, xl
  use M_fdata, only: nssh, rcutoff
  integer(c_long) ialp
  integer(c_long) iatom
  integer(c_long) ibeta
  integer(c_long) imu
  integer(c_long) in1, in2
  integer(c_long) ineigh
  integer(c_long) jatom
  integer(c_long) jbeta
  integer(c_long) jneigh
  integer(c_long) katom
  integer(c_long) kbeta
  integer(c_long) kneigh
  integer(c_long) num_neigh
  real(c_double) distance2
  real(c_double) range2
  real(c_double) rcutoff_i, rcutoff_j
  real(c_double), dimension (3) :: diff
  real(c_double), dimension (3) :: vec, vec1, vec2
  do ialp = 1,natoms
    num_neigh = 0
    do ineigh = 1, neighn(ialp)
      iatom = neigh_j(ineigh,ialp)
      ibeta = neigh_b(ineigh,ialp)
      if (.not. (iatom .eq. ialp .and. ibeta .eq. 0)) then
       in1 = imass(iatom)
       rcutoff_i = 0.0d0
       do imu = 1, nssh(in1)
         if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
       end do
       vec1(:) = ratom(:,iatom) + xl(:,ibeta)
       do jneigh = 1, neighn(ialp)
         jatom = neigh_j(jneigh,ialp)
         jbeta = neigh_b(jneigh,ialp)
         if (.not. (jatom .eq. ialp .and. jbeta .eq. 0) .and. (ineigh .lt. jneigh)) then      
           in2 = imass(jatom)
           rcutoff_j = 0.0d0
           do imu = 1, nssh(in2)
             if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
           end do
           vec2(:) = ratom(:,jatom) + xl(:,jbeta)
           distance2 = (vec2(1) - vec1(1))**2 + (vec2(2) - vec1(2))**2  + (vec2(3) - vec1(3))**2
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
           if (distance2 .le. range2) then
             num_neigh = num_neigh + 1
             neigh_comj(1,num_neigh,ialp) = iatom
             neigh_comb(1,num_neigh,ialp) = ibeta
             neigh_comj(2,num_neigh,ialp) = jatom
             neigh_comb(2,num_neigh,ialp) = jbeta
             diff = vec2 - vec1
             neigh_comm(num_neigh,ialp) = -9999
             do kneigh = 1, neighn(iatom)
               katom = neigh_j(kneigh,iatom)
               kbeta = neigh_b(kneigh,iatom)
               vec(:) = xl(:,kbeta) + ratom(:,katom) - ratom(:,iatom)
               if ((abs(vec(1) - diff(1)) .lt. 1.0d-4) .and. (abs(vec(2) - diff(2)) .lt. 1.0d-4) .and.  (abs(vec(3) - diff(3)) .lt. 1.0d-4)) neigh_comm(num_neigh,ialp) = kneigh
             end do
             if (neigh_comm(num_neigh,ialp) .eq. -9999) num_neigh = num_neigh - 1 
           end if
         end if
       end do
      end if
    end do
    neigh_comn(ialp) = num_neigh
  end do
  return
end subroutine common_neighbors
