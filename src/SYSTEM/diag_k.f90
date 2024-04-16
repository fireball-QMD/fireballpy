subroutine diag_k ( )
  use M_system
  implicit none
  integer ikpoint
  integer imu
  real, dimension (3) :: k_temp
  if (iqout .ne. 2) allocate (blowre (norbitals, norbitals, nkpoints))
  if (iqout .ne. 2 .and. icluster .ne. 1) allocate (blowim (norbitals, norbitals, nkpoints))
  allocate (bbnkre (norbitals, norbitals, nkpoints))
  if (icluster .eq. 0 .and. igamma .eq. 0) allocate (bbnkim (norbitals, norbitals, nkpoints))

  if (igamma .eq. 1) then
    k_temp(:) = special_k(:,1)
    call kspace_gamma (1, k_temp )

  else if (igamma .eq. 2) then
    k_temp(:) = special_k(:,1)
    call kspace_double_generalized(1, k_temp )

  else if (igamma .eq.3) then
    k_temp(:) = special_k(:,1)
    call kspace_double(1, k_temp )


  else
    do ikpoint = 1, nkpoints
      k_temp(:) = special_k(:,ikpoint)
      call kspace_double (ikpoint, k_temp )
    end do 
  end if


end subroutine diag_k

