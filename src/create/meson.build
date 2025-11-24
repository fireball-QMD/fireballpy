modules_create_srcs = [
  'coefficients.f90',
  'constants.f90',
  'dimensions.f90',
  'precision.f90',
  'quadrature.f90',
  'wavefunctions.f90',
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
]

integrals_srcs = [
  'append_string.f90',
  'exc1crho.f90',
  'gaunt.f90',
  'gauntReal.f90',
  'gleg.f90',
  'goverlap1c.f90',
  'mk_index.f90',
  'mk_indexDipX.f90',
  'mk_indexDipY.f90',
  'mk_indexPP.f90',
  'nuxc1crho.f90',
  'onecentervdip.f90',
  'onecenterxc.f90',
  'radial_integral.f90',
  'rescaled_psi.f90',
  'threecenter.f90',
  'threecenter_integral.f90',
  'twocenter.f90',
  'twocenter_integral.f90',
]

interpolaters_srcs = [
  'buildspline_1d.f90',
  'dpsiofr.f90',
  'interpolate2d.f90',
  'psiofr.f90',
  'rprimeofr.f90',
  'vnnaofr.f90',
  'vppofr.f90',
]

kinetic_srcs = [
  'bessel.f90',
  'kinetic.f90',
  'Pintegral.f90',
  'Tintegral.f90',
  'tree.f90',
]

readfiles_srcs = [
  'iofile.f90',
  'iofile2c.f90',
  'iofile2c_x.f90',
  'iofile2c_xcna.f90',
  'iofile3c.f90',
  'readcreate.f90',
  'readpsi.f90',
  'readtheory.f90',
  'readvnn.f90',
  'readvpp.f90',
]

xcc_srcs = [
  'bpgcxc.f90',
  'cepal.f90',
  'ceperley_alder.f90',
  'corlyp1c.f90',
  'corlyp2c.f90',
  'density_calc.f90',
  'dexc.f90',
  'dpotxc12.f90',
  'dpotxc12s.f90',
  'dvxc.f90',
  'dvxc3c.f90',
  'get_potxc1c.f90',
  'get_potxc2c.f90',
  'gga91_sr.f90',
  'ggacrad1c.f90',
  'ggacrad2c.f90',
  'ggaxrad1c.f90',
  'ggaxrad2c.f90',
  'ldaxc.f90',
  'lsdavwn.f90',
  'pbe.f90',
  'rho1c_store.f90',
  'rho2c_store.f90',
  'vxc.f90',
]

create_srcs = []
foreach src : modules_create_srcs
  create_srcs += 'MODULES' / src
endforeach
foreach src : xcc_srcs
  create_srcs += 'XC' / src
endforeach
foreach src : exx_srcs
  create_srcs += 'EXX' / src
endforeach
foreach src : integrals_srcs
  create_srcs += 'INTEGRALS' / src
endforeach
foreach src : interpolaters_srcs
  create_srcs += 'INTERPOLATERS' / src
endforeach
foreach src : kinetic_srcs
  create_srcs += 'KINETIC' / src
endforeach
foreach src : readfiles_srcs
  create_srcs += 'READFILES' / src
endforeach

# Deal with MPI
if mpi_dep.found()
  create_srcs += 'MPI.f90'
else
  create_srcs += 'noMPI.f90'
endif

create_srcs += 'create.f90'

# Build create executable
create_exe = executable('create.x',
  create_srcs,
  dependencies: [blas_dep, lapack_dep, mpi_dep],
  install: true,
  install_dir: py.get_install_dir() / 'fireballpy' / 'basis',
)
