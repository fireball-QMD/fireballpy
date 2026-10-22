! This file contains most of the parameters that you would want to
! change during general use.  These set maximum values for several
! things.  Larger values allow the program to handle bigger things,
! but also increase memory use and can cause slow performance if
! set way to big.
!
! Se eliminaron el 2026-10-22, por no usarlos nadie: norb_max,
! nrho_points, nz_points, wfmax_points, max_points_pp, max_points_na
! (solo dimensionaban los .inc muertos de funciones de onda,
! pseudopotenciales y potencial no neutro) y el flag 'superspline'.

! maximum number of species (different atomtypes)
        integer nspec_max
        parameter (nspec_max = 20)

! maximum number of shells 8 implies something like
! s, p, d, f, s*, p*, d*, and f* as one possibility
        integer nsh_max
        parameter (nsh_max = 8)

! maximum number of non-zero matrix orbitals for two-center interactions
        integer inter_max
        parameter (inter_max = 180)
