.. currentmodule:: fireballpy

.. _howtouse:

How to use :class:`Fireball`
============================

:class:`Fireball` as an :class:`ase.calculators.calculator.Calculator`
needs to be attached first to a valid :class:`ase.Atoms` object.
This may be done by setting the ``calc`` attribute:

    >>> from ase.build import molecule
    >>> from fireballpy import Fireball
    >>> atoms = molecule('CH4')
    >>> atoms.calc = Fireball(fdata='biology')

Do not worry about the ``fdata`` argument.
Its importance will become clear later.
We are now ready to perform some computations!
Let's get the energy, charges and forces:

    >>> atoms.get_potential_energy()
    -216.68505107072147
    >>> atoms.get_charges()
    array([-0.44507143,  0.11126785,  0.11126785,  0.11126786,  0.11126786])
    >>> atoms.get_forces()
    array([[ 1.73489876e-15,  5.49377383e-15, -1.99149954e-07],
    ...    [-2.81616200e-01, -2.81616200e-01, -2.81616191e-01],
    ...    [ 2.81616200e-01,  2.81616200e-01, -2.81616191e-01],
    ...    [-2.81616282e-01,  2.81616282e-01,  2.81616290e-01],
    ...    [ 2.81616282e-01, -2.81616282e-01,  2.81616290e-01]])

Nice! We got our first results for a methane molecule.

Now, if this is the first time you run this tutorial you may have noticed
what appeared to be a download before any computation took place.
This is because we requiere the precomputed integrals of our base functions.
How we select which basis we desire? Well that's were the ``fdata`` argument
comes into play!

In the next section we will discuss a bit more about *FData* and their crucial role.
