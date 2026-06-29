modules_wfs_srcs = [
  'begin_input.f90',
  'pp_storage.f90',
]

rcatms_srcs = [
  'clebsch_gordon.f90',
  'derf0.f90',
  'get_nlmesh_size.f90',
  'get_psi.f90',
  'get_vee.f90',
  'integrate_hpsi.f90',
  'pp.f90',
  'psirc.f90',
  'rcatms.f90',
  'rcatms_DMOL.f90',
  'rcatms_excite.f90',
  'rcatms_GSM.f90',
  'vshort.f90',
]

vnn_srcs = [
  'get_ewfmesh_size.f90',
  'get_wfmesh_size.f90',
  'psiofr.f90',
  'readpsi.f90',
  'readpsi_excite.f90',
  'vnn.f90',
  'vnn_excite.f90',
]

wfs_srcs = []
foreach src : modules_wfs_srcs
  wfs_srcs += 'MODULES' / src
endforeach
foreach src : rcatms_srcs
  wfs_srcs += 'RCATMS' / src
endforeach
foreach src : vnn_srcs
  wfs_srcs += 'VNN' / src
endforeach

# Build wfs library
wfs_lib = library(
  'fpy_wfs',
  wfs_srcs,
  dependencies: [libxc_dep, common_dep],
  install: true,
  install_dir: py.get_install_dir() / 'fireballpy.libs',
)
wfs_dep = declare_dependency(link_with: wfs_lib, sources: wfs_srcs)
