# Get included dirs
incdir_numpy = run_command(
  py,
  [
    '-c', 'import numpy; import os; os.chdir(os.path.join("..", "tools")); print(os.path.relpath(numpy.get_include()))',
  ],
  check: true,
).stdout().strip()
incdir_f2py = run_command(
  py,
  [
    '-c', 'import numpy.f2py; import os; os.chdir(os.path.join("..", "tools")); print(os.path.relpath(numpy.f2py.get_include()))',
  ],
  check: true,
).stdout().strip()
inc_numpy = include_directories(incdir_numpy)
inc_f2py = include_directories(incdir_f2py)

# Fortran object for f2py
fortranobject_c = incdir_f2py / 'fortranobject.c'
fortranobject_lib = static_library(
  '_fortranobject',
  fortranobject_c,
  dependencies: py_dep,
  include_directories: [inc_numpy, inc_f2py],
  gnu_symbol_visibility: 'hidden',
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_numpy, inc_f2py],
)

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
  fortranobject_c,
  fortran_args: fflags,
  dependencies: [fortranobject_dep, fireball_dep],
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

py.install_sources(
  [
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
  subdir: 'fireballpy',
)

subdir('basis')
subdir('bands')
subdir('outputs')
