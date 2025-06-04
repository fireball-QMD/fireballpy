# Build _fireball
# --quiet flag to supress stupid non-sense warnings
# one might have to take it out if it does not work to debug
fireball_source = custom_target(
  'fireball_module',
  input: 'fireball.f90', # No module hence no .f90 wrappers
  output: ['_fireballmodule.c', '_fireball-f2pywrappers.f'],
  command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fireball', '--lower', '--build-dir', '@OUTDIR@', '--quiet'],
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

# Create version.py
generate_version = custom_target(
  'generate-version',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: 'version.py',
  input: '../tools/gitversion.py',
  command: [
    '../tools/gitversion.py',
    '--meson-dist',
    '--write', 'fireballpy' / 'version.py',
  ],
  install_dir: py.get_install_dir() / 'fireballpy',
  install_tag: 'python-runtime',
)

py.install_sources([
    '__init__.py',
    '_fireball.pyi',
    '_errors.py',
    '_options.py',
    '_correction.py',
    'utils.py',
    'fdata.py',
    'atoms.py',
    'kpoints.py',
    'fireball.py',
    'orbitals.py',
    'ase.py',
    'fdata.toml',
  ],
  subdir: 'fireballpy'
)

subdir('basis')
subdir('bands')
subdir('outputs')
