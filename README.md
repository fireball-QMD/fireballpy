# FireballPy: minimal Fireball for Python

<!-- [![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://fireball-QMD.github.io/fireballpy/html/_static/lab) -->
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fireball-QMD/fireballpy/HEAD?labpath=examples/fireballpy_skeleton.ipynb) -->

- **Documentation:** <https://fireballpy.github.io>
- **Source code:** <https://github.com/fireball-QMD/fireballpy>
- **Fireball website:** <https://fireball-qmd.github.io>

Clone the repository to your computer:

    git clone https://github.com/fireball-QMD/fireballpy.git
    git clone git@github.com:fireball-QMD/fireballpy.git

## How to compile in local:

### UV

The first time you must run:

```bash
uv venv
uv pip install requirements/build_requirements.txt
uv sync
```

Changes should update automatically.

### Conda environment

Only the first time, create the environment:

```bash
conda create -f env.yml
conda activate fireballpy
```

Install the package:

```bash
python install.py --conda # For default gnu+OpenBLAS installation
python install.py --fast --conda # For non-float-respect gnu+OpenBLAS installation
python install.py --intel --conda # For ifx+MKL installation
python install.py --intel --fast --conda # For non-float-respect ifx+MKL installation
```

Note that OpenBLAS is already included.
The flag `--intel` may be substituted with `--intel-old` for use with old `ifort` and `icc` compilers.

### Virtual environment

Create the environment. Then install the dependencies:

```bash
python -m pip install -r requirements/build_requirements.txt
python -m pip install -r requirements/package_requirements.txt
```

Finally install the package:

```bash
python install.py # For default gnu+OpenBLAS installation
python install.py --fast # For non-float-respect gnu+OpenBLAS installation
python install.py --intel # For ifx+MKL installation
python install.py --intel --fast # For non-float-respect ifx+MKL installation
```

