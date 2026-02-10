import argparse
import os
import sys
import tempfile
from subprocess import Popen, PIPE, DEVNULL


def setup(meson, folder, args, env):
    p = Popen(meson + ['setup', folder] + args, stdout=PIPE, stderr=PIPE, env=env)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Setting up meson failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def comp(meson, folder, warn):
    p = Popen(meson + ['compile', '-C', folder], stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Compilation failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")
    elif warn:
        print(out.decode())


def install(meson, folder):
    p = Popen(meson + ['install', '-C', folder], stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if not (p.returncode == 0):
        raise RuntimeError(f"Installation failed!\n"
                           f"{out.decode()}\n"
                           f"{err.decode()}")


def version():
    p = Popen([sys.executable, os.path.join('tools', 'gitversion.py'),
               '--write', os.path.join('fireballpy', 'version.py')], stdout=PIPE, stderr=PIPE)
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
                        help="Use ifx, icx and MKL")
    parser.add_argument("--intel-old", action="store_true",
                        help="Use ifort, icc and MKL")
    parser.add_argument("--fast", action="store_true",
                        help="Apply non-float-respecting optimizations")
    parser.add_argument("--warn", action="store_true",
                        help="Print compiling info")
    args = parser.parse_args()

    # Check mpi
    ismpi = False
    try:
        p = Popen(['mpirun', '--version'], stdout=DEVNULL, stderr=DEVNULL)
        _ = p.communicate()
        ismpi = True
    except FileNotFoundError:
        ismpi = 0

    env = os.environ.copy()
    meson = [sys.executable, '-m', 'mesonbuild.mesonmain']
    setup_args = ['-Dpython.install_env=auto']
    if args.intel:
        env['CC'] = 'icx'
        env['FC'] = 'ifx'
        if ismpi:
            env['MPICC'] = "mpiicx"
            env['MPIFC'] = "mpiifx"
            setup_args += ['-Dmpi=true']
        setup_args += ['-Dblas=mkl-dynamic-ilp64-seq']
    elif args.intel_old:
        env['CC'] = 'icc'
        env['FC'] = 'ifort'
        if ismpi:
            env['MPICC'] = "mpicc"
            env['MPIFC'] = "mpifort"
            setup_args += ['-Dmpi=true']
        setup_args += ['-Dblas=mkl-dynamic-ilp64-seq']
    else:
        env['CC'] = 'gcc'
        env['FC'] = 'gfortran'
        if ismpi:
            env['CC'] = 'mpicc'
            env['FC'] = 'mpifort'
            setup_args += ['-Dmpi=true']
    if args.fast:
        setup_args += ['-Doptimization=3', '-Dbuildtype=custom']

    if args.build_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            setup(meson, tmp, setup_args, env)
            comp(meson, tmp, args.warn)
            install(meson, tmp)
    else:
        setup(meson, args.build_dir, setup_args, env)
        comp(meson, args.build_dir, args.warn)
        install(meson, args.build_dir)
    version()


if __name__ == '__main__':
    main()
