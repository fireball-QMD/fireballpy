module M_constants
  use, intrinsic :: iso_fortran_env, only: double => real64

  ! Bohr constant
  real(double), parameter :: abohr = 0.529177210544d0
  real(double), parameter :: eq2 = 14.39975d0
  real(double), parameter :: Hartree = 14.39975d0/abohr

  ! Rydberg/eV conversion
  real(double), parameter :: ryd = 13.6d0

  ! Boltzmann's constant
  real(double), parameter :: kb = 8.617343693082104d-5! 1.380662d-23 J/K * 6.24150974*10^18eV/J

  ! fovermp = f/mp = (1eV/1angstrom)/(mass proton) = angstrom/(fs**2).
  ! JPL 2000 I believe that this is incorrect after some lengthy discussion
  ! with Srini and Kurt.  The true conversion should be atomic mass units.
  !        real(double), parameter :: fovermp = 0.0095790d0
  real(double), parameter :: fovermp = 0.009648957597d0

  real(double), parameter :: kconvert = 11604.49558d0
  real(double), parameter :: pi = 3.141592653589793238462d0
  real(double), parameter :: spin = 2.0d0

  ! Gear algorithm constants_fireball
  !         integer, parameter :: gear_order = 5  !  (only use 2-7)
  !         real(double), dimension (0:gear_order) :: cfac
  ! JOM-info I change this so that I can run verlet (gear_order=2)
  integer, parameter :: gear_order = 2  !  (only use 2-7)
  real(double), dimension (0:5) :: cfac

  ! Kronecker delta
  real(double), dimension (3, 3) :: delk

  ! Levi-Civita parameter
  real(double), dimension (3, 3, 3) :: xlevi

  ! This is debatable, but they better be set high enough that the answers do not
  ! vary enough to concern you!  The forces are most susceptible to the norders.
  ! Higher order does not necessarily imply higher accuracy.  In fact too high
  ! of order can cause "ringing" (rapid occilations of the function).
  ! Negative numbers tell it to use a cubic spline of order abs(norder1)
  ! Splines provide better energy convservation and ovoid the ringing problem.
  ! Splines also avoid the problem with double numerical basis sets having
  ! unphysically low eigenvalues. 
  ! Note: what we used to call 6th order was real(double)ly 5th!
  ! We generally use -5 which puts 3 points on each side of where you are
  ! interpolating.  Total number of points=abs(norder1)+1.
  integer, parameter :: norder1 = -5
  ! What method do we use to interpolate 2-D grids
  integer, parameter :: D2intMeth = 1 ! 1 -> polynomials x then y
  ! All other methods sucked and were deleted
  ! twister data
  logical haveDorbitals

end module M_constants
