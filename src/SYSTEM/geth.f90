subroutine geth ()
  use iso_c_binding
  use M_system, only: natoms, imass, neighn, neigh_b, neigh_j, h_mat, s_mat, neighj_tot, neighb_tot, &
    & neighPP_j, neighPP_b, neighn, neighPPn, neighn_tot, numorb_max, vnl, hvec, svec, colvec, rowvec
  use M_fdata, only: num_orb
  implicit none
  integer(c_long) :: iatom, jatom, katom, imu, inu, in1, in2, mbeta, &
    & ineigh, numorb, jatom0, ineigh0, mbeta0, il, rowtemp, coltemp
  real(c_double), dimension (numorb_max, numorb_max) :: htemp, stemp

  numorb = 0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ineigh = 1, neighn_tot(iatom)
      jatom = neighj_tot(ineigh,iatom)
      in2 = imass(jatom)
      do imu = 1, num_orb(in1)
        do inu = 1, num_orb(in2)
          numorb = numorb + 1
        end do
      end do
    end do
  end do

  ! allocate linear vectors
  ! we want vectors with the column and row indices of the values in the
  ! linear vectors of the hamiltonian and overlap
  if (allocated(hvec)) deallocate(hvec)
  if (allocated(svec)) deallocate(svec)
  if (allocated(colvec)) deallocate(colvec)
  if (allocated(rowvec)) deallocate(rowvec)
  allocate(hvec(numorb), svec(numorb), colvec(numorb), rowvec(numorb))

  il = 1
  htemp = 0.0d0
  stemp = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)

    ! Row indices to be explored
    rowvec(il) = 0
    do katom = 1, iatom-1
      do imu = 1, num_orb(imass(katom))
        rowvec(il) = rowvec(il) + 1
      end do
    end do

    ! loop over total list of neighbors
    do ineigh = 1, neighn_tot(iatom)
      jatom = neighj_tot(ineigh,iatom)

      ! Col indices to be explored
      colvec(il) = 0
      do katom = 1, jatom-1
        do inu = 1, num_orb(imass(katom))
          colvec(il) = colvec(il) + 1
        end do
      end do

      mbeta = neighb_tot(ineigh,iatom)
      in2 = imass(jatom)
      htemp = 0.0d0
      stemp = 0.0d0

      ! loop over regular list of neighbors
      do ineigh0 = 1, neighn(iatom)
        jatom0 = neigh_j(ineigh0,iatom)
        mbeta0 = neigh_b(ineigh0,iatom)

        ! find identical neighbors
        if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              htemp(imu,inu) = htemp(imu,inu) + h_mat(imu,inu,ineigh0,iatom)
              stemp(imu,inu) = stemp(imu,inu) + s_mat(imu,inu,ineigh0,iatom)
            end do
          end do
        end if
      end do

      ! loop over PP list of neighbors
      do ineigh0 = 1, neighPPn(iatom)
        jatom0 = neighPP_j(ineigh0,iatom)
        mbeta0 = neighPP_b(ineigh0,iatom)

        ! find identical neighbors
        if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              htemp(imu,inu) = htemp(imu,inu) + vnl(imu,inu,ineigh0,iatom)
            end do
          end do
        end if
      end do

      ! write elements
      rowtemp = rowvec(il)
      do inu = 1, num_orb(in2)
        coltemp = colvec(il) + inu
        do imu = 1, num_orb(in1)
          colvec(il) = coltemp
          rowvec(il) = rowtemp + imu
          hvec(il) = htemp(imu,inu)
          svec(il) = stemp(imu,inu)
          il = il + 1
        end do
      end do
    end do
  end do
end subroutine geth
