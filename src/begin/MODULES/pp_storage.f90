        module pp_storage
         use precision
         
         integer npoints_short

         integer, dimension (:), allocatable :: npoints_nl

         real(kind=long) alpha
         real(kind=long) charge
         real(kind=long) drr_short
         real(kind=long) rrc_short
         real(kind=long) Zval

         real(kind=long), dimension (:), allocatable :: drr_nl
         real(kind=long), dimension (:), allocatable :: rrc_nl
         real(kind=long), dimension (:,:), allocatable :: r_nl
         real(kind=long), dimension (:), allocatable :: r_short
         real(kind=long), dimension (:,:), allocatable :: v_nl
         real(kind=long), dimension (:), allocatable :: v_short
        end module
