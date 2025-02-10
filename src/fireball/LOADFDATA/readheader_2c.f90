subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1, &
    & rc2, zmin, zmax, npseudo, cl_pseudo)
  use iso_c_binding
  implicit none
  integer(c_long), intent (in) :: interaction, iounit, nsh_max
  integer(c_long), intent (out) :: npseudo, numz
  real(c_double), intent (out) :: rc1, rc2, zmin, zmax
  real(c_double), intent (out), dimension (nsh_max) :: cl_pseudo
  integer(c_long) :: iline, issh, nucz1, nucz2

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
