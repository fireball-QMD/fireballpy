# FireballPy: minimal Fireball for Python

[![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://fireball-QMD.github.io/fireballpy/html/_static/lab)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fireball-QMD/fireballpy/HEAD?labpath=examples/fireballpy_skeleton.ipynb)

- **Documentation:** <https://fireball-QMD.github.io/fireballpy/html>, <https://fireballpy.github.io>
- **Source code:** <https://github.com/fireball-QMD/fireballpy>
- **Fireball website:** <https://fireball-qmd.github.io>

Clone the repository to your computer:

    git clone https://github.com/fireball-QMD/fireballpy
    git clone git@github.com:fireball-QMD/fireballpy
## How to compile in local:

First ensure you have the desired virtual Python environment activated (or conda environment).
Also ensure that you have either MKL or OpenBLAS installed in the system.
Then run:

```bash
python install.py --intel # For default ifx+MKL installation
python install.py --intel --fast # For non-float-respect ifx+MKL installation
python install.py # For default gnu+OpenBLAS installation
python install.py --fast # Fort non-float-respect gnu+OpenBLAS installation
```

To set a folder where to perform the build, set the `-b` (or `--build-dir`) flag.
Note this will, for the moment, only work on Linux
