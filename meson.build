project(
  'fireballpy',
  'fortran', 'c',
  version: '0.0.1',
  license: 'GPL-3',
  default_options: [
    'buildtype=debugoptimized',
    'fortran_std=legacy',
    'c_std=c17',
    'blas=openblas',
    'lapack=openblas',
  ],
)

# Check Fortran compiler and set flags
fc = meson.get_compiler('fortran')
fflags = ['-unroll', '-fPIC']
if fc.get_id() == 'gcc'
  fflags += ['-ffree-line-length-none', '-fdefault-real-8']
  if get_option('buildtype') == 'custom'
    fflags += ['-O3', '-ffast-math', '-fstack-arrays', '-fallow-store-data-races',
               '-fno-semantic-interposition', '-march=native']
  endif
elif fc.get_id() in ['intel', 'intel-llvm']
  fflags += ['-r8', '-fno-alias']
  if get_option('buildtype') == 'custom'
    fflags += ['-O3', '-fp-model=fast', '-xHost']
  else
    fflags += ['-fp-model=strict']
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

# Link LAPACK and BLAS (strat from scipy)
blas_name = get_option('blas')
lapack_name = get_option('lapack')
_args_blas_lapack = [] # Might be useful later

if blas_name == 'openblas'
  blas = dependency(['openblas', 'OpenBLAS'])
else # Different options
  blas = dependency(blas_name)
endif
blas_dep = declare_dependency(
  dependencies: blas,
  compile_args: _args_blas_lapack
)

if blas_name.startswith('mkl') or blas_name == 'openblas' # They bring LAPACK
  lapack = blas
else # Different options
  lapack = dependency(lapack_name)
endif
lapack_dep = declare_dependency(
  dependencies: lapack,
  compile_args: _args_blas_lapack
)

# Prepare f2py
incdir_numpy = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'],
                           check : true).stdout().strip()

incdir_f2py = run_command(py, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
                          check : true).stdout().strip()

fortranobject_lib = static_library('_fortranobject',
  incdir_f2py / 'fortranobject.c',
  dependencies : [py_dep],
  include_directories : [incdir_numpy, incdir_f2py]
)

fortranobject_dep = declare_dependency(
  link_with : fortranobject_lib,
  include_directories : [incdir_numpy, incdir_f2py]
)

# Build fireball
subdir('src')

# Build fireballpy
subdir('fireballpy')

# Test
test('molecule', py,
     args : [meson.current_source_dir() / 'test' / 'unittest_molecule.py'],
     env : env)
test('periodic', py,
     args : [meson.current_source_dir() / 'test' / 'unittest_periodic.py'],
     env : env,
     timeout: 0)
