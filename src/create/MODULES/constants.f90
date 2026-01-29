module constants
  use precision, only: wp
  implicit none
  public
  real(kind=wp), parameter :: abohr = 0.529177249_wp
  real(kind=wp), parameter :: abohr15 = 0.3849477153_wp
  real(kind=wp), parameter :: beta = 0.9_wp
  real(kind=wp), parameter :: eq2 = 14.39975_wp
  real(kind=wp), parameter :: Hartree = 14.39975_wp/abohr
  real(kind=wp), parameter :: pi = 3.141592653589793238462643_wp
  real(kind=wp), parameter :: ryd = 13.6057981_wp
  real(kind=wp), parameter :: tolerance = 1.0e-5_wp
  real(kind=wp), parameter :: inv4pi = 0.07957747154594767_wp
end module constants
