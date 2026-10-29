module M_constants
  use, intrinsic :: iso_fortran_env, only: dp => real64

  ! Bohr constant
  real(dp), parameter :: abohr = 0.529177210544_dp
  real(dp), parameter :: eq2 = 14.39975_dp
  real(dp), parameter :: Hartree = 14.39975_dp/abohr

  ! Rydberg/eV conversion
  real(dp), parameter :: ryd = 13.6_dp

  ! Boltzmann's constant
  real(dp), parameter :: kb = 8.617343693082104e-5_dp! 1.380662e-23_dp J/K * 6.24150974*10^18eV/J

  ! fovermp = f/mp = (1eV/1angstrom)/(mass proton) = angstrom/(fs**2).
  ! JPL 2000 I believe that this is incorrect after some lengthy discussion
  ! with Srini and Kurt.  The true conversion should be atomic mass units.
  !        real(dp), parameter :: fovermp = 0.0095790_dp
  real(dp), parameter :: fovermp = 0.009648957597_dp

  real(dp), parameter :: kconvert = 11604.49558_dp
  real(dp), parameter :: pi = 3.141592653589793238462_dp
  real(dp), parameter :: spin = 2.0_dp

  ! Gear algorithm constants_fireball
  !         integer, parameter :: gear_order = 5  !  (only use 2-7)
  !         real(dp), dimension (0:gear_order) :: cfac
  ! JOM-info I change this so that I can run verlet (gear_order=2)
  integer, parameter :: gear_order = 2  !  (only use 2-7)
  real(dp), dimension (0:5) :: cfac

  ! Kronecker delta
  real(dp), dimension (3, 3) :: delk

  ! Levi-Civita parameter
  real(dp), dimension (3, 3, 3) :: xlevi

  ! This is debatable, but they better be set high enough that the answers do not
  ! vary enough to concern you!  The forces are most susceptible to the norders.
  ! Higher order does not necessarily imply higher accuracy.  In fact too high
  ! of order can cause "ringing" (rapid occilations of the function).
  ! Negative numbers tell it to use a cubic spline of order abs(norder1)
  ! Splines provide better energy convservation and ovoid the ringing problem.
  ! Splines also avoid the problem with dp numerical basis sets having
  ! unphysically low eigenvalues. 
  ! Note: what we used to call 6th order was really 5th!
  ! We generally use -5 which puts 3 points on each side of where you are
  ! interpolating.  Total number of points=abs(norder1)+1.
  integer, parameter :: norder1 = -5
  ! What method do we use to interpolate 2-D grids
  integer, parameter :: D2intMeth = 1 ! 1 -> polynomials x then y
  ! All other methods sucked and were deleted
  ! twister data
  logical haveDorbitals

end module M_constants
