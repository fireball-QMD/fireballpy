        module pp_storage
         use, intrinsic :: iso_fortran_env, only: dp => real64
         
         integer npoints_short

         integer, dimension (:), allocatable :: npoints_nl

         real(kind=dp) alpha
         real(kind=dp) charge
         real(kind=dp) drr_short
         real(kind=dp) rrc_short
         real(kind=dp) Zval

         real(kind=dp), dimension (:), allocatable :: drr_nl
         real(kind=dp), dimension (:), allocatable :: rrc_nl
         real(kind=dp), dimension (:,:), allocatable :: r_nl
         real(kind=dp), dimension (:), allocatable :: r_short
         real(kind=dp), dimension (:,:), allocatable :: v_nl
         real(kind=dp), dimension (:), allocatable :: v_short
        end module
