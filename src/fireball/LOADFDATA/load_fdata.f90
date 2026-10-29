! Build the index tables and read the integrals.
!
! Everything that describes the FData -- how many species, which shells, which
! charges and radii -- is already in M_fdata, put there by 'set_fdata_info'
! from Python. Fortran no longer parses 'info.dat': the description travels as
! arguments, not as a file whose layout every reader had to count by hand.
subroutine load_fdata_tables()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nsh_max, nspecies, ind2c, icon3c, splineint_2c, &
    & isorpmax, isorpmax_xc, ME2c_max, initype, maxtype_2c
  implicit none
  integer :: in1, in2, in3, icount, isorp, interaction, interactions2c_max

  ! Compatibility
  isorpmax = nsh_max
  isorpmax_xc = nsh_max

  ! Allocate arrays
  call make_munu()
  call make_munuPP()
  call make_munuS()
  call make_munuDipY()
  call make_munuDipX()

  ! Load one centre
  call read_1c ()

  ! Prepare and load two centres
  icount = 0
  ind2c = 0
  interactions2c_max = (nsh_max+1)*3+8*nsh_max+9
  ! 'numz2c' and 'z2cmax' are gone: the z grid of every two centre table now
  ! travels inside its own math_interp_t, so there is nothing to book-keep.
  if (allocated(splineint_2c)) deallocate(splineint_2c)
  allocate (splineint_2c (ME2c_max, interactions2c_max, nspecies, nspecies))

  do interaction = 1, 20
    do isorp=initype(interaction),nsh_max*maxtype_2c(interaction)
      icount = icount + 1
      ind2c(interaction, isorp) = icount
    end do
    call read_2c(interaction)
  end do

  ! Prepare and load three centre
  if (allocated(icon3c)) deallocate(icon3c)
  allocate(icon3c(nspecies, nspecies, nspecies))
  do in1=1,nspecies
    do in2=1,nspecies
      do in3=1,nspecies
        icon3c(in1, in2, in3) = in3 + nspecies*(in2 - 1 + nspecies*(in1 - 1))
      end do
    end do
  end do
  do interaction = 1, 4
    if (interaction == 2) cycle
    call read_3c(interaction)
  end do
  call setterp_2d ()
end subroutine load_fdata_tables
