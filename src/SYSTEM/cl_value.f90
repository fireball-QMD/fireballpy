  ! This routine returns the Kleinman Bylander cl values for atom itype.
  ! The "raw" date is read in vnl.z1.z2.dat by  program read2c. We include up to
  ! 5 non-local (L values) of the pseudopotential. Usually, you will have 2
  ! (L = 0, 1) and sometimes 3 (L = 2).
subroutine cl_value (itype, cl)
  use M_system
  implicit none
  integer, intent (in) :: itype
  real, intent (out), dimension (numorb_max) :: cl
  integer imu
  integer index
  integer issh
  integer Lvalue
  integer Lmax
  cl(1:numorb_max) = 0.0d0
  ! We now loop though all shells, and create cl for each orbital.  For example,
  ! sp^3 has two shells; cl(1) = cl_PP(0) and cl(2) = cl(3) = cl(4) = cl_PP(1).
  index = 0
  do issh = 1, nsshPP(itype)
    Lvalue = lsshPP(issh,itype)
    Lmax = (2*Lvalue + 1)
    do imu = 1, Lmax
      index = index + 1
      cl(index) = cl_PP(issh,itype)
    end do
  end do
 
  ! Sanity check.
  if (index .ne. num_orbPP(itype)) then
    write (*,*) ' itype = ', itype
    write (*,*) ' index of orbitals for pseudopotential = ',index
    write (*,*) ' Program has num_orbPP = ', num_orbPP(itype)
    write (*,*) ' cl_value: index and num_orbPP DO NOT agree. '
    stop
  end if
  return
end
