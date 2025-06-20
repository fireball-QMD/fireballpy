subroutine diag_k ( )
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: iqout, igamma, Kscf, blowre, blowim, sm12_real, sm12_complex, special_k, norbitals, &
    & nkpoints, errno
  implicit none
  integer ikpoint
  real(double), dimension (3) :: k_temp

  if ((iqout .eq. 1) .or. (iqout .eq. 3)) then
    if (allocated(blowre)) deallocate(blowre)
    allocate (blowre (norbitals, norbitals, nkpoints))
    if (igamma .eq. 0) then
      if (allocated(blowim)) deallocate(blowim)
      allocate (blowim (norbitals, norbitals, nkpoints))
    end if
  end if


  if (igamma .eq. 0) then
    if (Kscf .eq. 1) then
      if(allocated(sm12_complex)) deallocate(sm12_complex)
      allocate (sm12_complex(norbitals,norbitals,nkpoints))
    end if
    do ikpoint = 1, nkpoints
      k_temp(:) = special_k(:,ikpoint)
      call kspace_double (ikpoint, k_temp )
      if (errno .ne. 0) return
    end do
  end if

  if (igamma .eq. 1) then
    if (Kscf .eq. 1) then
      if(allocated(sm12_real)) deallocate(sm12_real)
      allocate (sm12_real(norbitals,norbitals))
    end if
    call kspace_gamma ()
    if (errno .ne. 0) return
  end if

end subroutine diag_k
