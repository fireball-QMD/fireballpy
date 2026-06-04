! These common blocks store information about the wavefunctions
! You should not need to change any of this. But here it all is anyway.

        integer npoints (nsh_max, nspec_max)
        integer npoints_rho (nspec_max)
 
        integer nnrho
        integer nnz
 
        real(kind=wp) drho
        real(kind=wp) dz
        real(kind=wp) rhomax
        real(kind=wp) rhomin
        real(kind=wp) zmax
        real(kind=wp) zmin
 
        real(kind=wp) drr (nsh_max, nspec_max)               ! spacing between points
        real(kind=wp) drr_rho (nspec_max)                    ! spacing between points
        real(kind=wp) rr (wfmax_points, nsh_max, nspec_max)
        real(kind=wp) rr_rho (wfmax_points, nspec_max)
        real(kind=wp) rrc (nsh_max, nspec_max)               ! cutoff points
        real(kind=wp) rrc_rho (nspec_max)                    ! maximum cutoff points
        real(kind=wp) psi (wfmax_points, nsh_max, nspec_max)
        real(kind=wp) psi_spline (4,wfmax_points, nsh_max, nspec_max)
        real(kind=wp) rho2c (nrho_points, nz_points)
        real(kind=wp) rhop2c (nrho_points, nz_points)
        real(kind=wp) rhopp2c (nrho_points, nz_points)
        real(kind=wp) rhopz2c (nrho_points, nz_points)
        real(kind=wp) rhoz2c (nrho_points, nz_points)
        real(kind=wp) rhozz2c (nrho_points, nz_points)
 
        common /density_limits_rho/ drho, rhomin, rhomax, nnrho
        common /density_limits_z/ dz, zmin, zmax, nnz
        common /densities/ rho2c, rhop2c, rhopp2c, rhopz2c,rhoz2c, rhozz2c
        common /densparam/ drr_rho, rr_rho, rrc_rho, npoints_rho
        common /wavefunctions/ psi, psi_spline
        common /waveparam/ drr, rr, rrc, npoints
 
!$OMP   THREADPRIVATE(/density_limits_rho/,/density_limits_z/,/densities/)
 
