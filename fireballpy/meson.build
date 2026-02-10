# Build _fireball
# --quiet flag to supress stupid non-sense warnings
# one might have to take it out if it does not work to debug
fireball_source = custom_target(
  'fireball_module',
  input: 'fireball.f90', # No module hence no .f90 wrappers
  output: ['_fireballmodule.c', '_fireball-f2pywrappers.f'],
  command: [
    py,
    '-m', 'numpy.f2py',
    '@INPUT@',
    '-m', '_fireball',
    '--lower',
    '--build-dir', '@OUTDIR@',
    '--quiet',
  ],
)

py.extension_module(
  '_fireball',
  ['fireball.f90', fireball_source],
  fortran_args: fflags,
  dependencies: [fortranobject_dep, fireball_dep],
  link_language: 'fortran',
  install: true,
  subdir: 'fireballpy',
)

py.install_sources(
  [
    '__init__.py',
    '_fireball.pyi',
    '_errors.py',
    '_options.py',
    '_correction.py',
    'version.py',
    'utils.py',
    'fdata.py',
    'atoms.py',
    'kpoints.py',
    'fireball.py',
    'orbitals.py',
    'ase.py',
    'fdata.toml',
  ],
  subdir: 'fireballpy',
)

subdir('basis')
subdir('bands')
subdir('outputs')
