        module coefficients 
         use iso_fortran_env, only: dp => real64

! These are the factors which come from the coefficient in the Ylm
         real(kind=dp), dimension (:,:), allocatable :: clm

! These are the factors which come from the coefficients for the Legendre
! polynomials
!        real*8, dimension (:) :: ctheta
!        real*8, dimension (:) :: ctheta_weights

        end module
