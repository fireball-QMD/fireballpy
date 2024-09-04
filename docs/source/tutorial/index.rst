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

FireballPy consists of one class, :class:`Fireball` and one module
:mod:`retro`.

The :class:`Fireball` class provides an ASE compatible calculator to
perform DFT computations using ASE API.
This is possible thanks to the Fortran to C module conversion
performed by `f2py <https://numpy.org/doc/stable/f2py/>`_,
which although it provides with fantastic performance it does
also come with some shortcomings which we will explore later.

The module :mod:`retro` is dedicated to old time Fireball users
to ease with the adaptation of old projects and input files
to FireballPy.

.. toctree::
   :caption: Fireball calculator
   :maxdepth: 1

   howtouse

.. toctree::
   :caption: Coming from old Fireball
   :maxdepth: 1

   readingbas

.. toctree::
   :caption: Theory
   :maxdepth: 1

   mixer

.. toctree::
   :hidden:
   :caption: Extras

   ../license
