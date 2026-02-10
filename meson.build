project(
  'fireballpy',
  'fortran',
  'c',
  version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
  license: 'GPL-3.0-or-later',
  license_files: 'LICENSE',
  meson_version: '>=1.1.0',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c11',
    'fortran_std=f2008',
    'blas=openblas',
    'lapack=openblas',
    'mpi=false',
  ],
)

# Check Fortran compiler and set flags
fc = meson.get_compiler('fortran')
commonflags = ['-unroll', '-fPIC']
fflags = []
cflags = ['-Wno-misleading-indentation']

if fc.get_id() == 'gcc'
  # GNU
  fflags += ['-ffree-line-length-none']
  cflags += ['-Wno-uninitialized']
  if get_option('buildtype') == 'custom'
    commonflags += [
      '-ffast-math',
      '-fallow-store-data-races',
      '-fno-semantic-interposition',
      '-march=native',
    ]
    fflags += ['-fstack-arrays']
  endif
elif fc.get_id() in ['intel', 'intel-llvm']
  # Intel
  commonflags += ['-fno-alias']
  fflags += ['-free']
  if get_option('buildtype') == 'custom'
    commonflags += ['-fp-model=fast', '-xHost', '-diag-disable 10121']
  else
    commonflags += ['-fp-model=strict']
  endif
else
  error('No compiler found.')
endif

fflags += commonflags
cflags += commonflags

# Check python
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Set pkg-config environment
env = environment()
env.append('PATH', meson.current_source_dir())
env.append('PKG_CONFIG_PATH', get_option('libdir') / 'pkgconfig')

# Set compiler flags
add_global_arguments(fflags, language: 'fortran')
add_global_arguments(cflags, language: 'c')

# Find mpifort
if get_option('mpi')
  mpi_dep = dependency('mpi', language: 'fortran')
endif

# Link LAPACK and BLAS (strat from scipy)
blas_name = get_option('blas')
lapack_name = get_option('lapack')
_args_blas_lapack = [] # Might be useful later
if blas_name == 'openblas'
  # flexiblas is here for compatibility with fedora
  blas = dependency(['openblas', 'OpenBLAS', 'flexiblas'])
else
  blas = dependency(blas_name)
endif
blas_dep = declare_dependency(dependencies: blas, compile_args: _args_blas_lapack)

if blas_name.startswith('mkl') or blas_name == 'openblas'
  # They bring LAPACK
  lapack = blas
else
  lapack = dependency(lapack_name)
endif
lapack_dep = declare_dependency(dependencies: lapack, compile_args: _args_blas_lapack)

# Get included dirs
if py.full_path().contains(meson.current_source_dir())
  incdir_numpy = run_command(
    py,
    [
      '-c', 'import numpy; import os; print(os.path.relpath(numpy.get_include()))',
    ],
    check: true,
  ).stdout().strip()
  incdir_f2py = run_command(
    py,
    [
      '-c', 'import numpy.f2py; import os; print(os.path.relpath(numpy.f2py.get_include()))',
    ],
    check: true,
  ).stdout().strip()
else
  incdir_numpy = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()
  incdir_f2py = run_command(
    py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check: true,
  ).stdout().strip()
endif

fortranobject_lib = static_library(
  '_fortranobject',
  incdir_f2py / 'fortranobject.c',
  dependencies: [py_dep],
  include_directories: [incdir_numpy, incdir_f2py],
)

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [incdir_numpy, incdir_f2py],
)

# Build _fireball, _begin and create.x
subdir('src')

# Build fireballpy
subdir('fireballpy')

# Test
test(
  'molecule',
  py,
  args: [meson.current_source_dir() / 'test' / 'unittest_molecule.py'],
  env: env,
)
test(
  'periodic',
  py,
  args: [meson.current_source_dir() / 'test' / 'unittest_periodic.py'],
  env: env,
  timeout: 0,
)
