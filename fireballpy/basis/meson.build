# Build _wfs
# --quiet flag to supress stupid non-sense warnings
# one might have to take it out if it does not work to debug
wfs_source = custom_target(
  'wfs_module',
  input: 'wfs.f90', # No module hence no .f90 wrappers
  output: ['_wfsmodule.c', '_wfs-f2pywrappers.f'],
  command: [
    py,
    '-m', 'numpy.f2py',
    '@INPUT@',
    '-m', '_wfs',
    '--lower',
    '--build-dir', '@OUTDIR@',
    '--quiet',
  ],
)

py.extension_module(
  '_wfs',
  ['wfs.f90', wfs_source],
  fortran_args: fflags,
  dependencies: [fortranobject_dep, wfs_dep],
  install: true,
  subdir: 'fireballpy/basis',
)

py.install_sources(
  [
    '__init__.py',
    '__main__.py',
    '_wfs.pyi',
    'fdata.py',
  ],
  subdir: 'fireballpy/basis',
)
