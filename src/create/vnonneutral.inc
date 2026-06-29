! These common blocks store information about the (non)-neutral potential
! You should not need to change any of this. But here it all is anyway.

        integer npoints_na (0:nsh_max, nspec_max)
 
        real(kind=dp) drr_na (0:nsh_max, nspec_max)           ! spacing between points
        real(kind=dp) rr_na (max_points_na, 0:nsh_max, nspec_max)
        real(kind=dp) vnna (max_points_na, 0:nsh_max, nspec_max)
!        real(kind=dp) vnna_spline (4,max_points_na, 0:nsh_max, nspec_max)
 
        common /nonneutralpot/ rr_na, vnna!, vnna_spline
        common /nnaparam/ drr_na, npoints_na
