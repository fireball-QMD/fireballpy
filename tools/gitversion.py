#!/usr/bin/env python3
import os
import textwrap
import os.path
from subprocess import Popen, PIPE


# Reference: https://github.com/scipy/scipy/blob/main/tools/gitversion.py
def init_version():
    init = os.path.join(os.path.dirname(__file__), '../pyproject.toml')
    with open(init) as fid:
        data = fid.readlines()

    version_line = next(
        line for line in data if line.startswith('version =')
    )

    version = version_line.strip().split(' = ')[1]
    version = version.replace('"', '').replace("'", '')
    return version


def git_version(version):
    # Append last commit date and hash to dev version information, if available
    git_hash = ''
    try:
        # '%ad' with '--date=short' instead of '%aI': the latter needs git >= 2.2
        # and older gits (e.g. 1.8 on some clusters) emit the literal '%aI',
        # which used to end up verbatim in the version string.
        p = Popen(['git', 'log', '-1', '--date=short', '--format=%H %ad'], stdout=PIPE, stderr=PIPE,
                  cwd=os.path.dirname(__file__))
    except FileNotFoundError:
        pass
    else:
        out, err = p.communicate()
        if p.returncode == 0:
            fields = out.decode('utf-8').strip().replace('"', '').split()
            if len(fields) == 2:
                git_hash, git_date = fields[0], fields[1].replace('-', '')

                # Only attach git tag to development versions, and only if the
                # date really came out as a date: better no tag than a broken one
                if 'dev' in version and git_date.isdigit() and len(git_date) == 8:
                    version += f'+git{git_date}.{git_hash}'

    return version, git_hash


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--write', help="Save version to this file")
    parser.add_argument('--meson-dist', help='Output path is relative to MESON_DIST_ROOT', action='store_true')
    args = parser.parse_args()

    version, git_hash = git_version(init_version())

    template = textwrap.dedent(f'''
        """
        Module to expose more detailed version info for the installed `fireballpy`
        """
        version = "{version}"
        full_version = version
        short_version = version.split('.dev')[0]
        git_revision = "{git_hash}"
        release = 'dev' not in version and '+' not in version

        if not release:
            version = full_version
    ''')

    if args.write:
        outfile = args.write
        if args.meson_dist:
            outfile = os.path.join( os.environ.get('MESON_DIST_ROOT', ''), outfile)

        # Print human readable output path
        relpath = os.path.relpath(outfile)
        if relpath.startswith('.'):
            relpath = outfile

        with open(outfile, 'w') as f:
            f.write(template)
    else:
        print(version)
