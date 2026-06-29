        module x_exact
         use iso_fortran_env, only: dp => real64

         integer nnrhop

         real(kind=dp) drhop

         real(kind=dp), dimension (:), allocatable :: rpoint
         real(kind=dp), dimension (:,:,:), allocatable :: rprime
         real(kind=dp), dimension (:,:,:,:), allocatable :: rprime_spline
        end module
