! These common blocks store information about the wavefunctions
! You should not need to change any of this. But here it all is anyway.

        integer npoints (nsh_max, nspec_max)
        integer npoints_rho (nspec_max)
 
        integer nnrho
        integer nnz
 
        real(kind=long) drho
        real(kind=long) dz
        real(kind=long) rhomax
        real(kind=long) rhomin
        real(kind=long) zmax
        real(kind=long) zmin
 
        real(kind=long) drr (nsh_max, nspec_max)               ! spacing between points
        real(kind=long) drr_rho (nspec_max)                    ! spacing between points
        real(kind=long) rr (wfmax_points, nsh_max, nspec_max)
        real(kind=long) rr_rho (wfmax_points, nspec_max)
        real(kind=long) rrc (nsh_max, nspec_max)               ! cutoff points
        real(kind=long) rrc_rho (nspec_max)                    ! maximum cutoff points
        real(kind=long) psi (wfmax_points, nsh_max, nspec_max)
        real(kind=long) psi_spline (4,wfmax_points, nsh_max, nspec_max)
        real(kind=long) rho2c (nrho_points, nz_points)
        real(kind=long) rhop2c (nrho_points, nz_points)
        real(kind=long) rhopp2c (nrho_points, nz_points)
        real(kind=long) rhopz2c (nrho_points, nz_points)
        real(kind=long) rhoz2c (nrho_points, nz_points)
        real(kind=long) rhozz2c (nrho_points, nz_points)
 
        common /density_limits_rho/ drho, rhomin, rhomax, nnrho
        common /density_limits_z/ dz, zmin, zmax, nnz
        common /densities/ rho2c, rhop2c, rhopp2c, rhopz2c,rhoz2c, rhozz2c
        common /densparam/ drr_rho, rr_rho, rrc_rho, npoints_rho
        common /wavefunctions/ psi, psi_spline
        common /waveparam/ drr, rr, rrc, npoints
 
!$OMP   THREADPRIVATE(/density_limits_rho/,/density_limits_z/,/densities/)
 
