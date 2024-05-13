project(
  'fireballpy',
  'fortran', 'c',
  version: '0.0.1',
  license: 'GPL-3',
  default_options: [
    'buildtype=debugoptimized',
    'fortran_std=legacy',
    'c_std=c17',
    'blas=openblas'
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

# Set compiler flags
add_global_arguments(fflags, language: 'fortran')

# Link OpenBLAS or MKL
blas_name = get_option('blas')
if blas_name == 'openblas'
  blas_dep = dependency(['openblas', 'OpenBLAS'])
elif blas_name == 'mkl'
  blas_dep = dependency('mkl-dynamic-lp64-iomp')
endif

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

# Generate C module and Fortran wrappers
fireball_source = custom_target('fireball_module',
  input : ['src' / 'libf2py.f90'], # No module hence no .f90 wrappers
  output : ['_fireballmodule.c', '_fireball-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fireball', '--lower']
)

py.extension_module('_fireball',
  ['src' / 'libf2py.f90', fireball_source], # We use F90 so .f90 file is already included
  fortran_args : fflags,
  dependencies : [fortranobject_dep, blas_dep, fireball_dep],
  link_language : 'fortran',
  install : true,
  install_dir : meson.current_source_dir() / 'fireballpy' # Temporal
)
