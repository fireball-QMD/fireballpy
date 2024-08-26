subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1, &
    & rc2, zmin, zmax, npseudo, cl_pseudo)
  use M_constants, only: wp
  implicit none
  integer, intent (in) :: interaction, iounit, nsh_max
  integer, intent (out) :: npseudo, numz
  real(wp), intent (out) :: rc1, rc2, zmin, zmax
  real(wp), intent (out), dimension (nsh_max) :: cl_pseudo
  integer :: iline, issh, nucz1, nucz2

  do iline = 1, 9
      read (iounit,*)
  end do
  read (iounit,*) nucz1, rc1
  read (iounit,*) nucz2, rc2

  if (interaction .eq. 5) then
      read (iounit,*) npseudo
      read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
  end if

  read (iounit,*) zmax, numz
  zmin = 0.0d0
end subroutine readheader_2c
