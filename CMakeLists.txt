cmake_minimum_required(VERSION 3.0)

#set(ENV{FC} "ifort")

# Verificar si ifort está instalado
find_program(IFORT_EXECUTABLE ifort)

# Establecer la variable de entorno FC dependiendo del resultado
if(IFORT_EXECUTABLE)
    set(ENV{FC} "ifort")
    set(PATH_INCL ${MKLROOT}/include)
    set(PATH_MKL ${MKLROOT}/lib/intel64)
else()
    message(STATUS "ifort no encontrado. Usando gfortran.")
    find_program(GFORTRAN_EXECUTABLE gfortran)
    if(GFORTRAN_EXECUTABLE)
        set(ENV{FC} "gfortran")
    else()
        message(FATAL_ERROR "No se encontró ni ifort ni gfortran. Imposible compilar.")
    endif()
endif()


project(pyreball Fortran)

# Compilar programa en Fortran
add_executable(pyreball ${CMAKE_SOURCE_DIR}/src/pyreball.f90 
                        ${CMAKE_SOURCE_DIR}/src/diagonalize_matrix.f90
                        ${CMAKE_SOURCE_DIR}/src/variables.f90
                        )

if(IFORT_EXECUTABLE)
  #target_include_directories(pyreball PRIVATE ${MKL_INCL})
  set(FLAGS "-mkl")
  target_link_libraries(pyreball PRIVATE ${MKL_PATH}  ${FLAGS})
endif()
if(GFORTRAN_EXECUTABLE)
  set(FLAGS "-llapack -lblas")
  target_link_libraries(pyreball PRIVATE ${FLAGS})
endif()

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
