subroutine diag_k ( )
  use M_sytem
  implicit none
  integer ikpoint
  integer imu
  real, dimension (3) :: k_temp
  if (iqout .ne. 2) allocate (blowre (norbitals, norbitals, nkpoints))
  if (iqout .ne. 2 .and. icluster .ne. 1) allocate (blowim (norbitals, norbitals, nkpoints))
  allocate (bbnkre (norbitals, norbitals, nkpoints))
  if (icluster .ne. 1) allocate (bbnkim (norbitals, norbitals, nkpoints))
  allocate (eigen_k (norbitals, nkpoints))
  do ikpoint = 1, nkpoints
    k_temp(:) = special_k(:,ikpoint)
    !AQUI
    call kspace (Kscf, iqout, icluster,ikpoint, k_temp, nkpoints )
    if (ikpoint .eq. 1) then
      if (iwrteigen .eq. 1) write (19,*) nkpoints, norbitals_new
    end if
  end do ! do ikpoint
  return
end subroutine diag_k

