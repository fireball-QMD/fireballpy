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
    '_errors.py',
    '_options.py',
    '_correction.py',
    'utils.py',
    'fdata.py',
    'atoms.py',
    'kpoints.py',
    'fireball.py',
    'ase.py',
    'fdata.toml',
  ],
  subdir: 'fireballpy'
)

subdir('bands')
subdir('outputs')
