project(
  'fireballpy',
  'fortran',
  'c',
  version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
  license: 'GPL-3',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c11',
    'fortran_std=f2008',
    'blas=openblas',
    'lapack=openblas',
  ],
)

# Check Fortran compiler and set flags
fc = meson.get_compiler('fortran')
fflags = ['-unroll', '-fPIC']
cflags = ['-unroll', '-fPIC', '-Wno-misleading-indentation']
if fc.get_id() == 'gcc'
  fflags += ['-ffree-line-length-none', '-Wno-maybe-uninitialized']
  if get_option('buildtype') == 'custom'
    # unitialized warning could be problematic as gnu is different from intel
    fflags += [
      '-ffast-math',
      '-fstack-arrays',
      '-fallow-store-data-races',
      '-fno-semantic-interposition',
      '-march=native',
      '-Wno-uninitialized',
    ]
    cflags += [
      '-ffast-math',
      '-fallow-store-data-races',
      '-fno-semantic-interposition',
      '-march=native',
      '-Wno-uninitialized',
    ]
  endif
elif fc.get_id() in ['intel', 'intel-llvm']
  fflags += ['-fno-alias', '-free']
  cflags += ['-fno-alias']
  if get_option('buildtype') == 'custom'
    fflags += ['-fp-model=fast', '-xHost', '-diag-disable 10121']
    cflags += ['-fp-model=fast', '-xHost', '-diag-disable 10121']
  else
    fflags += ['-fp-model=strict']
    cflags += ['-fp-model=strict']
  endif
else
  error('No compiler found.')
endif

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
#mpi_dep = dependency('mpi', language: 'fortran', required: false)

# Link LAPACK and BLAS (strat from scipy)
blas_name = get_option('blas')
lapack_name = get_option('lapack')
_args_blas_lapack = [] # Might be useful later

if blas_name == 'openblas'
  blas = dependency(['openblas', 'OpenBLAS', 'flexiblas'])
else
  # Different options
  blas = dependency(blas_name)
endif
blas_dep = declare_dependency(dependencies: blas, compile_args: _args_blas_lapack)

if blas_name.startswith('mkl') or blas_name == 'openblas'
  # They bring LAPACK
  lapack = blas
else
  # Different options
  lapack = dependency(lapack_name)
endif
lapack_dep = declare_dependency(dependencies: lapack, compile_args: _args_blas_lapack)

# Prepare f2py
incdir_numpy = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()

incdir_f2py = run_command(
  py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true,
).stdout().strip()

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
