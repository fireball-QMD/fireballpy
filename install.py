import argparse
import os
import stat
import sys
import tempfile
from subprocess import Popen, PIPE


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
    parser.add_argument("--mpi", action="store_true",
                        help="Compile with mpi")
    parser.add_argument("--warn", action="store_true",
                        help="Print compiling info")
    parser.add_argument("--conda", action="store_true",
                        help="Make links to work with conda envs")
    args = parser.parse_args()

    env = os.environ.copy()
    meson = [sys.executable, '-m', 'mesonbuild.mesonmain']
    setup_args = ['-Dpython.install_env=auto']
    if args.intel:
        env['CC'] = 'icx'
        env['FC'] = 'ifx'
        if args.mpi:
            env['MPICC'] = "mpiicx"
            env['MPIFC'] = "mpiifx"
        setup_args += ['-Dblas=mkl-dynamic-ilp64-seq']
    elif args.intel_old:
        env['CC'] = 'icc'
        env['FC'] = 'ifort'
        if args.mpi:
            env['MPICC'] = "mpicc"
            env['MPIFC'] = "mpifort"
        setup_args += ['-Dblas=mkl-dynamic-ilp64-seq']
    else:
        env['CC'] = 'gcc'
        env['FC'] = 'gfortran'
        if args.mpi:
            env['CC'] = 'mpicc'
            env['FC'] = 'mpifort'
    if args.fast:
        setup_args += ['-Doptimization=3', '-Dbuildtype=custom']
    if args.mpi:
        setup_args += ['-Dmpi=true']

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

    if args.conda:
        pyver = '.'.join(sys.version.split('.')[0:2])
        lib_folder = os.path.join(sys.prefix, sys.platlibdir)
        fpy_folder = os.path.join(lib_folder, 'python' + pyver, 'site-packages', 'fireballpy.libs')
        for lib in ['libfireball.so', 'libbegin.so']:
            if os.path.isfile(os.path.join(lib_folder, lib)):
                os.remove(os.path.join(lib_folder, lib))
            os.symlink(os.path.join(fpy_folder, lib), os.path.join(lib_folder, lib))

        bin_folder = os.path.join(sys.prefix, 'bin')
        if os.path.isfile(os.path.join(bin_folder, 'fdata')):
            os.remove(os.path.join(bin_folder, 'fdata'))
        with open(os.path.join(bin_folder, 'fdata'), 'x') as file:
            file.write('#!' + sys.executable + '\n')
            with open(os.path.join(lib_folder, 'python' + pyver, 'site-packages', 'fireballpy', 'basis', 'fdata.py'), 'r') as fp:
                fdata = fp.read()
            file.write(fdata)
        os.chmod(os.path.join(bin_folder, 'fdata'), stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)


if __name__ == '__main__':
    main()
