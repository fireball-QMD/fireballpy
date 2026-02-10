# Build _begin
# --quiet flag to supress stupid non-sense warnings
# one might have to take it out if it does not work to debug
begin_source = custom_target(
  'begin_module',
  input: 'begin.f90', # No module hence no .f90 wrappers
  output: ['_beginmodule.c', '_begin-f2pywrappers.f'],
  command: [
    py,
    '-m', 'numpy.f2py',
    '@INPUT@',
    '-m', '_begin',
    '--lower',
    '--build-dir', '@OUTDIR@',
    '--quiet',
  ],
)

py.extension_module(
  '_begin',
  ['begin.f90', begin_source],
  fortran_args: fflags,
  dependencies: [fortranobject_dep, begin_dep],
  link_language: 'fortran',
  install: true,
  subdir: 'fireballpy/basis',
)

py.install_sources(
  [
    '__init__.py',
    '__main__.py',
    '_begin.pyi',
    'fdata.py',
  ],
  subdir: 'fireballpy/basis',
)
