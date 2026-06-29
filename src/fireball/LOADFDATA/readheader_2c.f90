subroutine readheader_2c (interaction, iounit, in1, in2, nsh_max, numz, rc1, &
    & rc2, zmin, zmax, npseudo, cl_pseudo)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nssh, Qref
  implicit none
  integer, intent (in) :: interaction, iounit, nsh_max, in1, in2
  integer, intent (out) :: npseudo, numz
  real(double), intent (out) :: rc1, rc2, zmin, zmax
  real(double), intent (out), dimension (nsh_max) :: cl_pseudo
  integer :: iline, issh, nucz1, nucz2

  ! TODO: for all
  if (interaction == 13) then
    do iline = 1, 9
        read (iounit,*)
    end do
    read (iounit,*) nucz1, rc1
    read (iounit,*) nucz2, rc2
    read (iounit,*) zmax, numz
    zmin = 0.0d0
    return
  end if

  read (iounit,*) nucz1, nucz2
  read (iounit,*) rc1, rc2
  read (iounit,*) zmin, zmax, numz
  if (interaction == 5) then
    read (iounit,*) npseudo
    read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
  else if (interaction >= 6 .and. interaction <= 8) then
    read (iounit,*)
    read (iounit,*) (Qref(issh, in1), issh = 1, nssh(in1))
    read (iounit,*) (Qref(issh, in2), issh = 1, nssh(in2))
  end if
  read (iounit,*)
end subroutine readheader_2c
