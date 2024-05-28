subroutine diag_k ( )
  use M_system
  implicit none
  integer ikpoint
  integer imu
  real*8, dimension (3) :: k_temp
  if (allocated(bbnkre)) deallocate(bbnkre)
  allocate (bbnkre (norbitals, norbitals, nkpoints))

  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    if (allocated(blowre)) deallocate(blowre)
    allocate (blowre (norbitals, norbitals, nkpoints))
    if (igamma .eq. 0) then
      if (allocated(blowim)) deallocate(blowim)
      allocate (blowim (norbitals, norbitals, nkpoints))
    end if
  end if

  if (icluster .eq. 0 .and. igamma .eq. 0) then
      if (allocated(bbnkim)) deallocate(bbnkim)
      allocate (bbnkim (norbitals, norbitals, nkpoints))
  end if
 
  if (igamma .eq. 0) then
    if (Kscf .eq. 1) then
      if(allocated(sm12_complex)) deallocate(sm12_complex)
      allocate (sm12_complex(norbitals,norbitals,nkpoints))
    end if
    do ikpoint = 1, nkpoints      
      k_temp(:) = special_k(:,ikpoint)
     call kspace_double (ikpoint, k_temp )
    end do 
  end if
 
  if (igamma .eq. 1) then  
    if (Kscf .eq. 1) then
       if(allocated(sm12_real)) deallocate(sm12_real)
       allocate (sm12_real(norbitals,norbitals))
    end if
    k_temp(:) = special_k(:,1)
    call kspace_gamma (1, k_temp )
  end if

  if (igamma .eq. 2) then
    k_temp(:) = special_k(:,1)
    call kspace_double_generalized(1, k_temp )
  end if

  if (igamma .eq. 3) then
    k_temp(:) = special_k(:,1)
    call kspace_double(1, k_temp )
  end if


end subroutine diag_k

