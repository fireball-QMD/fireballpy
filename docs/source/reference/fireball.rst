.. currentmodule:: fireballpy

.. _fireballpy-api-fireball:

Fireball low-level API
======================

.. autoclass:: AtomSystem
   :members: set_coords, set_cell, update_coords
   :no-inherited-members:

.. autoclass:: FDataFiles
   :members: load_fdata, get_charges_method, get_correction
   :no-inherited-members:

.. autoclass:: KPoints
   :members: reduce_kpts, set_kpoints
   :no-inherited-members:

.. autoclass:: BaseFireball
   :members: compute_charges, compute_energy, compute_eigenvalues, compute_forces, update_coords
   :no-inherited-members:
