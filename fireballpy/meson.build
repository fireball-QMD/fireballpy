# Create version.py
generate_version = custom_target(
  'generate-version',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: 'version.py',
  input: '../tools/gitversion.py',
  command: ['../tools/gitversion.py', '--meson-dist', '--write', 'fireballpy' / 'version.py'],
  install_dir: py.get_install_dir() / 'fireballpy',
  install_tag: 'python-runtime',
)

python_sources = [
  '__init__.py',
  'utils.py',
  'fdata.py',
  'infodat.py',
  'retro.py',
  'fireball.py',
]

py.install_sources(
  python_sources,
  subdir: 'fireballpy'
)
