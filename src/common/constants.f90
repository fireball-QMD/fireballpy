module constants
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
private
integer, parameter, public :: ERR_LEN = 256
integer, parameter, public :: FNAME_LEN = 255
integer, parameter, public :: PATH_LEN = 4096
real(kind=dp), parameter, public :: ABOHR = 0.5291772105638411_dp
real(kind=dp), parameter, public :: ABOHR1_5 = 0.38494767331059465_dp
real(kind=dp), parameter, public :: ABOHR3 = 0.14818471118724036_dp
real(kind=dp), parameter, public :: ABOHR4 = 0.07841597211427229_dp
real(kind=dp), parameter, public :: ABOHR5 = 0.04149594538708256_dp
real(kind=dp), parameter, public :: ABOHR6 = 0.02195870862964584_dp
real(kind=dp), parameter, public :: ABOHR8 = 0.006149064682626328_dp
real(kind=dp), parameter, public :: ABOHR13 = 0.0002551612522519003_dp
real(kind=dp), parameter, public :: INVABOHR = 1.8897261258369282_dp
real(kind=dp), parameter, public :: BETA = 0.9_dp
real(kind=dp), parameter, public :: EQ2 = 14.399645351950548_dp
real(kind=dp), parameter, public :: HARTREE = 27.211386024367243_dp
real(kind=dp), parameter, public :: RYD = 13.605693012183622_dp
real(kind=dp), parameter, public :: TOLERANCE = 1.0e-8_dp
real(kind=dp), parameter, public :: PI = 3.141592653589793_dp
real(kind=dp), parameter, public :: TWOPI = 6.283185307179586_dp
real(kind=dp), parameter, public :: INV4PI = 0.07957747154594767_dp
real(kind=dp), parameter, public :: SQINV4PI = 0.28209479177387814_dp
real(kind=dp), parameter, public :: INVSQ2 = 0.7071067811865475_dp

! Begin specific
real(kind=dp), parameter, public :: TOLERANCE_SCF = 1.0e-5_dp
real(kind=dp), parameter, public :: A0_EFF(0:3) = [2.0_dp, 1.0_dp, 0.8_dp, 0.7_dp]

! Names
character(1), parameter, public :: LNAMES(0:3) = ["s", "p", "d", "f"]
character(32), parameter, public :: SNAMES(0:0) = ["s"]
character(32), parameter, public :: PNAMES(-1:1) = ["p_{y}", "p_{z}", "p_{x}"]
character(32), parameter, public :: DNAMES(-2:2) = ["d_{xy}     ", "d_{yz}     ", "d_{z^2}    ", "d_{xz}     ", "d_{x^2-y^2}"]
character(32), parameter, public :: FNAMES(-3:3) = ["f_{y(3x^2-y^2)}", "f_{xyz}        ", "f_{yz^2}       ", &
  &                                                 "f_{z^3}        ", &
  &                                                 "f_{xz^2}       ", "f_{z(x^2-y^2)} ", "f_{x(x^2-3y^2)}"]
end module constants
