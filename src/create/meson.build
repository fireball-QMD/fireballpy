basic_srcs = [
  'wavefunctions.f90',
  'potentials.f90',
  'pseudopotentials.f90'
]

modules_create_srcs = [
  'coefficients.f90',
  'dimensions.f90',
  'quadrature.f90',
  'x_exact.f90',
]

exx_srcs = [
  'clebsch_gordon.f90',
  'factorial.f90',
  'setup_clm.f90',
  'x_1c.f90',
  'xatom_2c.f90',
  'xatom_2c_integral.f90',
  'xatom_2c_rprime.f90',
  'xontopl_2c.f90',
  'xontopl_2c_integral.f90',
  'xontopl_2c_rprime.f90',
  'xontopr_2c.f90',
  'xontopr_2c_integral.f90',
  'xontopr_2c_rprime.f90',
  'rprimeofr.f90',
]

integrals_srcs = [
  'gaunt.f90',
  'gauntReal.f90',
  'gleg.f90',
  'mk_index.f90',
  'mk_indexPP.f90',
  'rescaled_psi.f90',
  'onecentervdip.f90',
  'onecenter.f90',
  'radial_integral.f90',
  'threecenter.f90',
  'threecenter_integral.f90',
  'twocenter.f90',
]

kinetic_srcs = [
  'bessel.f90',
  'kinetic.f90',
  'Pintegral.f90',
  'Tintegral.f90',
  'tree.f90',
]

readfiles_srcs = [
  'iofile2c.f90',
  'iofile2c_x.f90',
  'iofile3c.f90',
  'readcreate.f90',
  'readpsi.f90',
  'readtheory.f90',
  'readvnn.f90',
  'readvpp.f90',
]

create_srcs = ['create.f90']
foreach src : basic_srcs
  create_srcs += src
endforeach
foreach src : modules_create_srcs
  create_srcs += 'MODULES' / src
endforeach
foreach src : exx_srcs
  create_srcs += 'EXX' / src
endforeach
foreach src : integrals_srcs
  create_srcs += 'INTEGRALS' / src
endforeach
foreach src : kinetic_srcs
  create_srcs += 'KINETIC' / src
endforeach
foreach src : readfiles_srcs
  create_srcs += 'READFILES' / src
endforeach

# Deal with MPI
if get_option('mpi')
  create_srcs += 'MPI.f90'
  create_exe = executable(
    'create.x',
    create_srcs,
    dependencies: [blas_dep, lapack_dep, mpi_dep, libxc_dep, common_dep],
    install: true,
    install_dir: py.get_install_dir() / 'fireballpy' / 'basis',
  )
else
  create_srcs += 'noMPI.f90'
  create_exe = executable(
    'create.x',
    create_srcs,
    dependencies: [blas_dep, lapack_dep, libxc_dep, common_dep],
    install: true,
    install_dir: py.get_install_dir() / 'fireballpy' / 'basis',
  )
endif
