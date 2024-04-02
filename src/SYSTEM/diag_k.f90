subroutine diag_k ( )
  use M_system
  implicit none
  integer ikpoint
  integer imu
  real, dimension (3) :: k_temp
  if (iqout .ne. 2) allocate (blowre (norbitals, norbitals, nkpoints))
  if (iqout .ne. 2 .and. icluster .ne. 1) allocate (blowim (norbitals, norbitals, nkpoints))
  allocate (bbnkre (norbitals, norbitals, nkpoints))
  if (icluster .ne. 1) allocate (bbnkim (norbitals, norbitals, nkpoints))

  if (gamma .eq. 1) then
    k_temp(:) = special_k(:,1)
    call kspace_gamma (1, k_temp )
  else
    do ikpoint = 1, nkpoints
      k_temp(:) = special_k(:,ikpoint)
      call kspace_double (ikpoint, k_temp )
    end do 
  end if
end subroutine diag_k

