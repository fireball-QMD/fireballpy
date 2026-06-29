module constants
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  integer, parameter, public :: err_len = 256
  integer, parameter, public :: fname_len = 255
  integer, parameter, public :: path_len = 4096
  real(kind=dp), parameter, public :: abohr = 0.5291772105638411_dp
  real(kind=dp), parameter, public :: abohr1_5 = 0.38494767331059465_dp
  real(kind=dp), parameter, public :: abohr3 = 0.14818471118724036_dp
  real(kind=dp), parameter, public :: abohr4 = 0.07841597211427229_dp
  real(kind=dp), parameter, public :: abohr5 = 0.04149594538708256_dp
  real(kind=dp), parameter, public :: abohr8 = 0.006149064682626328_dp
  real(kind=dp), parameter, public :: abohr13 = 0.0002551612522519003_dp
  real(kind=dp), parameter, public :: invabohr = 1.8897261258369282_dp
  real(kind=dp), parameter, public :: beta = 0.9_dp
  real(kind=dp), parameter, public :: eq2 = 14.399645351950548_dp
  real(kind=dp), parameter, public :: hartree = 27.211386024367243_dp
  real(kind=dp), parameter, public :: ryd = 13.605693012183622_dp
  real(kind=dp), parameter, public :: tolerance = 1.0e-5_dp
  real(kind=dp), parameter, public :: pi = 3.141592653589793_dp
  real(kind=dp), parameter, public :: twopi = 6.283185307179586_dp
  real(kind=dp), parameter, public :: inv4pi = 0.07957747154594767_dp
  real(kind=dp), parameter, public :: sqinv4pi = 0.28209479177387814_dp
  real(kind=dp), parameter, public :: invsq2 = 0.7071067811865475_dp

  ! Names
  character(32), parameter, public :: snames(0:0) = ["s"]
  character(32), parameter, public :: pnames(-1:1) = ["p_{y}", "p_{z}", "p_{x}"]
  character(32), parameter, public :: dnames(-2:2) = ["d_{xy}     ", "d_{yz}     ", "d_{z^2}    ", "d_{xz}     ", "d_{x^2-y^2}"]
  character(32), parameter, public :: fnames(-3:3) = ["f_{y(3x^2-y^2)}", "f_{xyz}        ", "f_{yz^2}       ", &
    &                                                 "f_{z^3}        ", &
    &                                                 "f_{xz^2}       ", "f_{z(x^2-y^2)} ", "f_{x(x^2-3y^2)}"]
end module constants
