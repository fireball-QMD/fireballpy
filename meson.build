project('fireballpy', 'fortran', 'c',
        version: '0.0.1',
        license: 'GPL-3',
        default_options: ['optimization=3'])

# Check Fortran compiler
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
    flags = ['-fdefault-real-8', '-unroll', '-fPIC', '-ffree-line-length-512', '-fimplicit-none']
else
    error('No compiler found.')
endif

# Set compiler flags
add_global_arguments(flags, language: 'fortran')

# Link OpenBLAS
blas = dependency(['openblas', 'OpenBLAS'])
blas_dep = declare_dependency(dependencies: blas)

# Build fireball
subdir('src')
