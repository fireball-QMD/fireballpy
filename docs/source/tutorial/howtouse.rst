.. currentmodule:: fireballpy

.. _howtouse:

How to use :class:`Fireball`
============================

:class:`Fireball` as an :class:`ase.calculators.calculator.Calculator`
needs to be attached first to a valid :class:`ase.Atoms` object.
This may be done by setting the ``calc`` attribute:

    >>> from ase import Atoms
    >>> from fireballpy import Fireball
    >>> atoms = Atoms('CH4',
    ...               positions=[( 0.0000,  0.0000,  0.0000),
    ...                          ( 0.6276,  0.6276,  0.6276),
    ...                          ( 0.6276, -0.6276, -0.6276),
    ...                          (-0.6276,  0.6276, -0.6276),
    ...                          (-0.6276, -0.6276,  0.6276)])
    >>> atoms.calc = Fireball(fdata='biology')

Do not worry about the ``fdata`` argument.
Its importance will become clear later.
We are now ready to perform some computations!
Let's get the energy, charges and forces:

    >>> atoms.get_potential_energy()
    -208.34774557250086
    >>> atoms.get_charges()
    array([-0.4479097, 0.11197742, 0.11197743, 0.11197743, 0.11197742])
    >>> atoms.get_forces()
    array([[ 4.30489932e-15, -2.74825372e-15, -2.00305221e-07],
    ...    [-2.22048842e-01, -2.22048842e-01, -2.22048833e-01],
    ...    [-2.22048924e-01,  2.22048924e-01,  2.22048933e-01],
    ...    [ 2.22048924e-01, -2.22048924e-01,  2.22048933e-01],
    ...    [ 2.22048842e-01,  2.22048842e-01, -2.22048833e-01]])

Nice! We got our first results for `this methane molecule <https://cccbdb.nist.gov/exp2x.asp?casno=74828&charge=0>`_.

Now, if this is the first time you run this tutorial you may have noticed
what appeared to be a download before any computation took place.
This is because we requiere the precomputed integrals of our base functions.
How we select which basis we desire? Well that's were the ``fdata`` argument
comes into play!

In the next section we will discuss a bit more about *FData* and their crucial role.
