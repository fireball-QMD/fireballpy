cmake_minimum_required(VERSION 3.0)

set(ENV{FC} "ifort")

project(pyreball Fortran)

set(BLA_VENDOR Intel10_64lp)
#set(MKLROOT "/opt/intel/oneapi/mkl/latest")
set(PATH_INCL ${MKLROOT}/include)
set(PATH_MKL ${MKLROOT}/lib/intel64)
set(MKL_FLAGS "-mkl")

# Compilar programa en Fortran
add_executable(pyreball ${CMAKE_SOURCE_DIR}/src/pyreball.f90 
                        ${CMAKE_SOURCE_DIR}/src/diagonalize_matrix.f90
                        ${CMAKE_SOURCE_DIR}/src/variables.f90
                        )

#target_include_directories(pyreball PRIVATE ${MKL_INCL})
target_link_libraries(pyreball PRIVATE ${MKL_PATH}  ${MKL_FLAGS})


#add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/src/fortran_suma.so
#                   COMMAND f2py -c ${CMAKE_SOURCE_DIR}/src/suma.f90
#                   ${CMAKE_SOURCE_DIR}/src/variables.f90 
#                   -m ${CMAKE_SOURCE_DIR}/src/f2py.f90 --f90flags='-fPIC'
#                   COMMAND mv ${CMAKE_SOURCE_DIR}/src/f2py ${CMAKE_SOURCE_DIR}/bin/
#                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
#                   COMMENT "Compilando la subrutina Fortran con f2py"
#                   VERBATIM)
#
#add_custom_target(fortran_suma ALL DEPENDS ${CMAKE_SOURCE_DIR}/src/fortran_suma.so)
