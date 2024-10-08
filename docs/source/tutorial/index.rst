.. currentmodule:: fireballpy

.. _user_guide:

*********************
FireballPy User Guide
*********************

**Release:** |release|
**Date:** |today|

FireballPy is an ASE calculator which implements the Fireball code
to perform DFT computations.
As we will see during this tutorial, speed is a key feature of Fireball,
which is achieved by the use of precomputed integral tables.

In this user guide we will provide a gentle introduction with examples
to the various options Fireball gives us.
A more detailed documentation may be found in :ref:`fireballpy-api`.
Please note that, although we will explain part of the ASE API for this
tutorial, this should not be taken as a replacement for
`ASE reference guide <https://wiki.fysik.dtu.dk/ase/>`_ which we
deeply recommend reading.

Structure of the package
========================

The :class:`Fireball` class provides an ASE compatible calculator to
perform DFT computations using ASE API.
This is possible thanks to the Fortran wrappers generated
by `f2py <https://numpy.org/doc/stable/f2py/>`_,
which although it provides with fantastic performance it does
also come with some shortcomings which we will explore later.

.. toctree::
   :caption: Fireball calculator
   :maxdepth: 1

   howtouse
   fdata

.. toctree::
   :caption: Theory
   :maxdepth: 1

   mixer
   charges

.. toctree::
   :hidden:
   :caption: Extras

   ../license
