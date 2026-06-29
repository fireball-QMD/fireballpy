common_srcs = [
  'constants.f90',
  'indices.f90',
  'math.f90',
  'xc.f90',
  'utils.f90'
]

common_lib = library(
  'fpy_common',
  common_srcs,
  dependencies: [libxc_dep],
  install: true,
  install_dir: py.get_install_dir() / 'fireballpy.libs',
)
common_dep = declare_dependency(link_with: common_lib, sources: common_srcs)
