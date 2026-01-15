        module constants
         use precision, only: wp
         
         real(kind=wp), parameter :: abohr = 0.529177249d0 
         real(kind=wp), parameter :: abohr15 = 0.3849477153d0
         real(kind=wp), parameter :: beta = 0.9d0
         real(kind=wp), parameter :: eq2 = 14.39975d0
         real(kind=wp), parameter :: Hartree = 14.39975d0/abohr 
         real(kind=wp), parameter :: pi = 3.141592653589793238462643d0
         real(kind=wp), parameter :: ryd = 13.6057981d0 
         real(kind=wp), parameter :: tolerance = 1.0d-5
        end module
