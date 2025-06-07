import argparse
import os
import sys
import subprocess
import tempfile


def setup(meson, folder, args, env):
    p = subprocess.Popen(meson + ['setup', folder] + args,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         env=env)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Setting up meson failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def comp(meson, folder):
    p = subprocess.Popen(meson + ['compile', '-C', folder],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    #print(out.decode())
    if not (p.returncode == 0):
        raise RuntimeError(f"Compilation failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def install(meson, folder):
    p = subprocess.Popen(meson + ['install', '-C', folder],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Installation failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def version():
    p = subprocess.Popen([sys.executable,
                          os.path.join('tools', 'gitversion.py'),
                          '--write',
                          os.path.join('fireballpy', 'version.py')],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Creating a local version.py file failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--build-dir", type=str,
                        help="Custom build folder (default temp directory)")
    parser.add_argument("--intel", action="store_true",
                        help="Use ifx and MKL")
    parser.add_argument("--fast", action="store_true",
                        help="Apply non-float-respecting optmisations")
    args = parser.parse_args()

    # Check mpi
    ismpi = 0
    #try:
    #    p = subprocess.Popen(['mpirun', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #    _ = p.communicate()
    #    ismpi = True
    #except FileNotFoundError:
    #    ismpi = 0

    env = os.environ.copy()
    meson = [sys.executable, '-m', 'mesonbuild.mesonmain']
    setup_args = ['-Dpython.install_env=auto']
    if args.intel:
        if ismpi:
            env['CC'] = f"mpiicx"
            env['FC'] = f"mpiifx"
        else:
            env['CC'] = 'icx'
            env['FC'] = 'ifx'
        setup_args += ['-Dblas=mkl-dynamic-lp64-iomp']
    else:
        if ismpi:
            env['CC'] = f"mpicc"
            env['FC'] = f"mpifort"
        else:
            env['CC'] = 'gcc'
            env['FC'] = 'gfortran'
    if args.fast:
        setup_args += ['-Doptimization=3', '-Dbuildtype=custom']

    if args.build_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            setup(meson, tmp, setup_args, env)
            comp(meson, tmp)
            install(meson, tmp)
    else:
        setup(meson, args.build_dir, setup_args, env)
        comp(meson, args.build_dir)
        install(meson, args.build_dir)
    version()


if __name__ == '__main__':
    main()
