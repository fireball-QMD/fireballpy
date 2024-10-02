.. currentmodule:: fireballpy

.. _howtouse:

How to use :class:`Fireball`
============================

:class:`Fireball` as an :class:`ase.calculators.calculator.Calculator`
needs to be attached first to a valid :class:`ase.Atoms` object.
This may be done by setting the ``calc`` attribute:

    >>> from ase.builder import molecule
    >>> from fireballpy import Fireball
    >>> atoms = molecule('CH4')
    >>> atoms = Atoms('CH4')
    >>> atoms.calc = Fireball(fdata='biology')

Do not worry about the ``fdata`` argument.
Its importance will become clear later.
We are now ready to perform some computations!
Let's get the energy, charges and forces:

    >>> atoms.get_potential_energy()
    -208.56448015852973
    >>> atoms.get_charges()
    array([-0.44507143,  0.11126785,  0.11126785,  0.11126786,  0.11126786])
    >>> atoms.get_forces()
    array([[ 3.18824910e-14,  3.49360951e-14, -4.16910533e-07],
    ...    [-1.56317853e+00, -1.56317853e+00, -1.56317854e+00],
    ...    [ 1.56317853e+00,  1.56317853e+00, -1.56317854e+00],
    ...    [-1.56317875e+00,  1.56317875e+00,  1.56317875e+00],
    ...    [ 1.56317875e+00, -1.56317875e+00,  1.56317875e+00]])

Nice! We got our first results for a methane molecule.

Now, if this is the first time you run this tutorial you may have noticed
what appeared to be a download before any computation took place.
This is because we requiere the precomputed integrals of our base functions.
How we select which basis we desire? Well that's were the ``fdata`` argument
comes into play!

In the next section we will discuss a bit more about *FData* and their crucial role.
