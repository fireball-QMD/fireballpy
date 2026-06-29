! These common blocks store information about the (non)-neutral potential
! You should not need to change any of this. But here it all is anyway.

        integer npoints_pp (0:nsh_max, nspec_max)
 
        real(kind=dp) cl_pp (nsh_max, nspec_max)
        real(kind=dp) drr_pp (nsh_max, nspec_max)           ! spacing between points
        real(kind=dp) rr_pp (max_points_pp, nsh_max, nspec_max)
        real(kind=dp) rrc_pp (nsh_max, nspec_max)
        real(kind=dp) vpp (max_points_pp, nsh_max, nspec_max)

!        real(kind=dp) vpp_spline (4,max_points_pp, nsh_max, nspec_max)
        real(kind=dp) rcPP (nspec_max)                                                                                             
        common /pp/ rr_pp, vpp!, vpp_spline
        common /pp_param/ cl_pp, drr_pp, rrc_pp, npoints_pp, rcPP

