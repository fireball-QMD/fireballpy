module constants
  use precision, only: wp
  implicit none
  public
  real(kind=wp), parameter :: abohr = 0.5291772105638411_wp
  real(kind=wp), parameter :: abohr15 = 0.38494767331059465_wp
  real(kind=wp), parameter :: abohr3 = 0.14818471118724036_wp
  real(kind=wp), parameter :: abohr4 = 0.07841597211427229_wp
  real(kind=wp), parameter :: abohr5 = 0.04149594538708256_wp
  real(kind=wp), parameter :: invabohr = 1.8897261258369282_wp
  real(kind=wp), parameter :: beta = 0.9_wp
  real(kind=wp), parameter :: eq2 = 14.399645351950548_wp
  real(kind=wp), parameter :: hartree = 27.211386024367243_wp
  real(kind=wp), parameter :: pi = 3.141592653589793_wp
  real(kind=wp), parameter :: ryd = 13.605693012183622_wp
  real(kind=wp), parameter :: tolerance = 1.0e-5_wp
  real(kind=wp), parameter :: inv4pi = 0.07957747154594767_wp
end module constants
