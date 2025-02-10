        module coefficients 
         use precision

! These are the factors which come from the coefficient in the Ylm
         real(kind=long), dimension (:,:), allocatable :: clm

! These are the factors which come from the coefficients for the Legendre
! polynomials
!        real*8, dimension (:) :: ctheta
!        real*8, dimension (:) :: ctheta_weights

        end module
