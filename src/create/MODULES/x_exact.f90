        module x_exact
         use precision, only: wp

         integer nnrhop

         real(kind=wp) drhop

         real(kind=wp), dimension (:), allocatable :: rpoint
         real(kind=wp), dimension (:,:,:), allocatable :: rprime
         real(kind=wp), dimension (:,:,:,:), allocatable :: rprime_spline
        end module
