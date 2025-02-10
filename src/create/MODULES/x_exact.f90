        module x_exact
         use precision

         integer nnrhop

         real(kind=long) drhop

         real(kind=long), dimension (:), allocatable :: rpoint
         real(kind=long), dimension (:,:,:), allocatable :: rprime
         real(kind=long), dimension (:,:,:,:), allocatable :: rprime_spline
        end module
