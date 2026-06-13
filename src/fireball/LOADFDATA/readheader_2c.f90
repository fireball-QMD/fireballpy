subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1, &
    & rc2, zmin, zmax, npseudo, cl_pseudo)
  use, intrinsic :: iso_fortran_env, only: double => real64
  implicit none
  integer, intent (in) :: interaction, iounit, nsh_max
  integer, intent (out) :: npseudo, numz
  real(double), intent (out) :: rc1, rc2, zmin, zmax
  real(double), intent (out), dimension (nsh_max) :: cl_pseudo
  integer :: iline, issh, nucz1, nucz2, nskip

  ! The new exchange-correlation 2-center files (vxc_2c, vxc_2c_ol, vxc_2c_or =
  ! interactions 6,7,8) are written with one fewer header line than the standard
  ! 2-center format (they lack the leading separator line), so we must skip 8
  ! header lines for them instead of 9. Skipping 9 misaligns the reader and the
  ! zmax/numz + spline data of the XC-ontop integral get read from the wrong lines.
  nskip = 9
  if (interaction .ge. 6 .and. interaction .le. 8) nskip = 8
  do iline = 1, nskip
      read (iounit,*)
  end do
  read (iounit,*) nucz1, rc1
  read (iounit,*) nucz2, rc2

  ! Only the non-local pseudopotential (vnl, interaction 5) carries the cl
  ! coefficients. interaction 6 (vxc_2c) does NOT, so it must not be read here.
  if (interaction .eq. 5) then
      read (iounit,*) npseudo
      read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
  end if

  read (iounit,*) zmax, numz
  zmin = 0.0d0
end subroutine readheader_2c
