! This file contains most of the parameters that you would want to
! change during general use.  These set maximum values for several
! things.  Larger values allow the program to handle bigger things,
! but also increase memory use and can cause slow performance if 
! set way to big.

! maximum number of species (different atomtypes)
        integer nspec_max
        parameter (nspec_max = 20)
        integer norb_max
 
! maximum number of shells 8 implies something like
! s, p, d, f, s*, p*, d*, and f* as one possibility
        integer nsh_max
        parameter (nsh_max = 8)
 
! maximum number of non-zero matrix orbitals for two-center interactions
!       integer inter2c_max
!       parameter (inter2c_max = 24)

! maximum number of non-zero matrix orbitals for two-center interactions
!       integer inter3c_max
!       parameter (inter3c_max = 40)

! maximum number of non-zero matrix orbitals for two-center interactions
        integer inter_max
!        parameter (inter_max = 97)
        parameter (inter_max = 180)

! This parameter gives the maximum number of radial points contained in the
! wavefunction array.
        integer nrho_points
        parameter (nrho_points = 1601)

        integer nz_points
        parameter (nz_points = 2501)

        integer wfmax_points
        parameter (wfmax_points = 1501)

! This parameter gives the maximum number of radial points contained in the
! (non)-neutral potential array.
        integer max_points_pp
        parameter (max_points_pp = 1501)

! This parameter gives the maximum number of radial points contained in the
! (non)-neutral potential array.
        integer max_points_na
        parameter (max_points_na = 1501)

! Superspline anyone?
        logical superspline
        parameter (superspline = .true.)


