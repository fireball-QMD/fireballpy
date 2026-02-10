
"""
Module to expose more detailed version info for the installed `fireballpy`
"""
version = "0.7.1.dev0+git20260112.167d421059c4822640e24cfdd82ffbec8e2d762b"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "167d421059c4822640e24cfdd82ffbec8e2d762b"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
