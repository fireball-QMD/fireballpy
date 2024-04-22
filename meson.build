project('fireballpy', 'fortran', 'c',
        version: '0.0.1',
        license: 'GPL-3')

# Check Fortran compiler
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
    flags = ['-unroll', '-fPIC', '-ffree-line-length-512', '-fimplicit-none']
else
    error('No compiler found.')
endif

# Set compiler flags
add_global_arguments(flags, language: 'fortran')

# Link OpenBLAS
blas_dep = dependency(['openblas', 'OpenBLAS'])

# Build fireball
subdir('src')
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

fireball_source = custom_target('fireball_module',
  input : ['src' / 'libf2py.f90'],
  output : ['_fireballmodule.c', '_fireball-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fireball', '--lower']
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('_fireball',
  ['src' / 'libf2py.f90', fireball_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : [py_dep, fireball_dep],
  install : true,
  install_dir : meson.current_source_dir() / 'fireballpy'
)
