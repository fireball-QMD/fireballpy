subroutine readheader_2c (interaction, iounit, in1, in2, nsh_max, numz, rc1, &
    & rc2, zmin, zmax, npseudo, cl_pseudo)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nssh, Qref, TWOCENTER_VNL, TWOCENTER_KINETIC, TWOCENTER_VXC_0, TWOCENTER_VXC_L, TWOCENTER_VXC_R
  implicit none
  integer, intent (in) :: interaction, iounit, nsh_max, in1, in2
  integer, intent (out) :: npseudo, numz
  real(double), intent (out) :: rc1, rc2, zmin, zmax
  real(double), intent (out) :: cl_pseudo(nsh_max)
  integer :: iline, issh, nucz1, nucz2, tempn1, tempn2
  real(double), allocatable :: tempq1(:), tempq2(:)

  ! TODO: for all
  if (interaction == TWOCENTER_KINETIC) then
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
  if (interaction == TWOCENTER_VNL) then
    read (iounit,*) npseudo
    read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
  else if (interaction == TWOCENTER_VXC_0 .or. interaction == TWOCENTER_VXC_L .or. interaction == TWOCENTER_VXC_R) then
    ! We use this to check everything is consistent with 1c
    read (iounit,*) tempn1, tempn2
    if (tempn1 /= nssh(in1) .or. tempn2 /= nssh(in2)) error stop
    allocate(tempq1(tempn1), tempq2(tempn2))
    read (iounit,*) (tempq1(issh), issh = 1, tempn1)
    if (any(tempq1 /= Qref(:,in1))) error stop
    read (iounit,*) (tempq2(issh), issh = 1, tempn2)
    if (any(tempq2 /= Qref(:,in2))) error stop
  end if
  read (iounit,*)
end subroutine readheader_2c
