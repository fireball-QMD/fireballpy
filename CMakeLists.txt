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
add_executable( pyreball 
 

  ${CMAKE_SOURCE_DIR}/src/M_constants.f90
  ${CMAKE_SOURCE_DIR}/src/M_fdata.f90 
  ${CMAKE_SOURCE_DIR}/src/M_system.f90
 
  
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/make_munuDipX.f90 
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/make_munu.f90 
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/make_munuS.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/make_munuDipY.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/make_munuPP.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/read_1c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/readheader_2c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/readdata_2c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/read_2c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/read_3c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/readdata_3c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/readheader_3c.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/setterp_2d.f90
  ${CMAKE_SOURCE_DIR}/src/LOADFDATA/load_fdata.f90
 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/load_system.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/allocate_system.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/anderson.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/broyden.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/louie.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/pulay.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mixer.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/fermie.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/denmat.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/build_rho.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/kspace_withnok.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/diag_k.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/diag_error.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/scf_loop.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/getenergy.f90                  
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/get_ewald.f90                  
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_usr.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/interpolate_1d.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recoverC.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_mcweda.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_olsxc_1c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/unocentros.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_olsxc_on.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_olsxc_off.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/average_ca_rho.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_ca_2c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_ca_3c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_lr.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/buildh.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_ca_3c_dip.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_lr_dip.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_3c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_ca_2c_dip.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_sVNL.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_2c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_2c_PP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_1c_vdip.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_3c_PP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/initneighbors.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/num_neigh_tot.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/backnay.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/neighbors_pairs.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/common_neighbors.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/common_neighborsPP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_olsxc_on.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/build_ca_olsxc_on.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_olsxc_off.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/epsilon.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/deps2cent.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentros.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosS.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/build_olsxc_off.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/average_ca_rho.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/trescentros.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/trescentrosS.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/smoother.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosPP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosDipY.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosDipX.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/cl_value.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mpairnay.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/rotated.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/rotate.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/common_neighbors.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_S.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_2c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/interpolate_2d.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_3c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosPP.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_PP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/rotatePP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/twister.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/twisterd.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/chooser.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/chooserd.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/rotatedPP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/makeDmatPP.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/makeDmat.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_2cDipX.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/recover_2cDipY.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosDipX.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/doscentrosDipY.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/cepal.f90
  
 # ${CMAKE_SOURCE_DIR}/src/SYSTEM/scf_loop.f90
             
  ${CMAKE_SOURCE_DIR}/src/pyreball.f90
  )
  #${CMAKE_SOURCE_DIR}/src/diagonalize_matrix.f90
  #${CMAKE_SOURCE_DIR}/src/variables.f90
  #)

if(IFORT_EXECUTABLE)
  set(FLAGS "-mkl")
  target_link_libraries(pyreball PRIVATE ${MKL_PATH}  ${FLAGS})
endif()
#if(GFORTRAN_EXECUTABLE)
#  set(FLAGS "-llapack -lblas")
#  target_link_libraries(pyreball PRIVATE ${FLAGS})
#endif()
#
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
