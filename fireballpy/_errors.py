from __future__ import annotations
from numpy.linalg import LinAlgError


def type_check(obj, objtype, name, msg=''):
    if not isinstance(obj, objtype):
        raise TypeError(f"Parameter ``{name}`` must be of type ``{objtype}``{msg}.")


def raise_fb_error(fb_errno):
    # Linear Algebra errors
    if fb_errno < 0:
        raise LinAlgError("Diagonalization not successful: "
                          f"{-fb_errno} off-diagonal elements failed "
                          "to converge. "
                          "Perhaps reducing 'beta' and 'mix_order' in "
                          "'mixer_kws' would help")

    # Denmat error
    if fb_errno == 1:
        raise RuntimeError("Total charge is not being conserved, most "
                           "probably for rounding errors. "
                           "Try running again")
