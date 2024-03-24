cmake_minimum_required(VERSION 3.0)

# Verificar si ifort está instalado
find_program(IFORT_EXECUTABLE ifort)

# Establecer la variable de entorno FC dependiendo del resultado
if(IFORT_EXECUTABLE)
  set(FFLAGS "-O3 -r8  -fno-alias -unroll -fPIC")
  set(ENV{FC} "ifort ${FFLAGS}")
  set(MKL_FLAGS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm")
  set(MKLROOT "/opt/intel/oneapi/mkl/latest")
  set(MKL_PATH "${MKLROOT}/lib/intel64/")
else()
  message(STATUS "ifort no encontrado. Usando gfortran.")
  find_program(GFORTRAN_EXECUTABLE gfortran)
  if(GFORTRAN_EXECUTABLE)
    set(ENV{FC} "gfortran")
  else()
    message(FATAL_ERROR "No se encontró ni ifort ni gfortran. Imposible compilar.")
  endif()
endif()

project(fireballpy Fortran)

add_library(fireballpy
 
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
 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/kspace_gamma.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/kspace_double.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/getforces.f90  
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mpairnay.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_2c.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_2c_PP.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_olsxc_on.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_olsxc_2c.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_2c.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_2c_dip.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_3c.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_3c_PP.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_3c.f90    
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_3c_dip.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_lr.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_lr_dip.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dassemble_ca_olsxc_3c.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/assemble_F.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/Dtrescentros.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/DtrescentrosS.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/deps3center.f90   
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/initamat.f90 
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/dipole_proyection.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mulliken_charges.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mulliken_dipole_charges.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/lowdin_charges.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/neighbors.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/initboxes.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/allocate_system.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/anderson.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/broyden.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/louie.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/pulay.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/mixer.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/fermie.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/denmat.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/build_rho.f90
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
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/neighbors.f90
  ${CMAKE_SOURCE_DIR}/src/SYSTEM/neighborsPP.f90
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
  )

if(IFORT_EXECUTABLE)
  target_link_libraries(fireballpy PUBLIC ${MKL_PATH}  ${MKL_FLAGS})
endif()
