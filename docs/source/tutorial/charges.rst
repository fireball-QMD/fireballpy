.. _charges_methods:

**********************
Charges mixing methods
**********************

Options are:

+--------------------+------------------------------+---------------------------+
| Method             | Value                        | Notes                     |
+====================+==============================+===========================+
| Automatic          | None                         | FData dependent           |
| (default)          |                              |                           |
+--------------------+------------------------------+---------------------------+
| Mulliken           | 'mulliken'                   | [Mulliken]_               |
+--------------------+------------------------------+---------------------------+
| Mulliken Dipole    | 'mulliken_dipole'            | [MullikenDip]_            |
+--------------------+------------------------------+---------------------------+
| Mulliken Dipole    | 'mulliken_dipole_preserving' |                           |
| Preserving         |                              |                           |
+--------------------+------------------------------+---------------------------+
| Lowdin             | 'lowdin'                     | [Lowdin]_                 |
+--------------------+------------------------------+---------------------------+
| Weighted Lowdin    | 'weighted_lowdin'            | [WeightLowdin]_. Reduces  |
|                    |                              | to Lowdin if the basis    |
|                    |                              | has no empty shells.      |
+--------------------+------------------------------+---------------------------+
| Stationary charges | 'stationary_charges'         | Charges that make the     |
|                    |                              | energy stationary,        |
|                    |                              | ∂E/∂Q = µ. See below.     |
+--------------------+------------------------------+---------------------------+

Stationary charges and ``fix_shells``
-------------------------------------

With ``charges_method='stationary_charges'`` the shell charges are obtained by
solving a linear system that makes the energy stationary, so that the
Hellmann-Feynman-like forces are consistent with the energy surface.

Empty shells (zero neutral charge: polarization d shells, excited s'/p' of
double bases) create near-null directions of that system and must be kept
fixed at their neutral charge. This is controlled with the ``fix_shells``
parameter of the calculator (ignored by all other charge methods):

+---------------------+---------------------------------------------------------------+
| ``fix_shells``      | Effect                                                        |
+=====================+===============================================================+
| ``'auto'``/``None`` | Fix every shell with zero neutral charge (default).           |
+---------------------+---------------------------------------------------------------+
| ``'d'``             | Fix only l=2 shells.                                          |
+---------------------+---------------------------------------------------------------+
| ``'none'``          | All shells free.                                              |
+---------------------+---------------------------------------------------------------+
| 0/1 mask            | Explicit per-shell control, one entry per shell in global     |
|                     | shell order (shells of atom 1, then atom 2, ...); 1 = fixed.  |
+---------------------+---------------------------------------------------------------+

The mask has one entry per shell, atom by atom in the order of the ``Atoms``
object, and within each atom its shells in the order defined by the FData.
For example, for a water molecule ordered ``OHH`` with a double basis,
O(s, p, s', p') + H(s, s') + H(s, s'):

.. code-block:: python

   from fireballpy import Fireball

   # default, equivalent to fix_shells='auto'
   atoms.calc = Fireball(fdata='custom', fdata_path='...',
                         charges_method='stationary_charges')

   # explicit mask: fix the excited shells (same as 'auto' here)
   #                O: s  p  s' p'  H: s  s'  H: s  s'
   atoms.calc = Fireball(fdata='custom', fdata_path='...',
                         charges_method='stationary_charges',
                         fix_shells=[0, 0, 1, 1, 0, 1, 0, 1])

A mask with the wrong length (it must match the total number of shells of
the system) or with values other than 0/1 raises a ``ValueError``.

.. warning::
   With a basis that has no empty shells (all s, p occupied) the default
   ``'auto'`` fixes nothing and the stationary system may converge to
   unphysical charges in some geometries. In that case consider providing
   an explicit mask.

Note that this is different from the ``fix_charges`` parameter, which
freezes *all* shell charges and skips the SCF loop entirely.

References
----------

.. [Mulliken] Reference for Mulliken.
              Must add.

.. [MullikenDip] Reference for Mulliken Dipole.
                 Must add.

.. [Lowdin] Reference for Lowdin.
            Must add.

.. [WeightLowdin] Reference for NPA.
         Must add.

