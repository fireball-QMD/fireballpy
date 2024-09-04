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
                      positions=[( 0.0000,  0.0000,  0.0000),
                                 ( 0.6276,  0.6276,  0.6276),
                                 ( 0.6276, -0.6276, -0.6276),
                                 (-0.6276,  0.6276, -0.6276),
                                 (-0.6276, -0.6276,  0.6276)])
    >>> atoms.calc = Fireball()

We are now ready to perform some computations!
Let's get the energy and the forces:

    >>> atoms.get_potential_energy()
    -216.27239440823337
    >>> atoms.get_forces()
    array([[ 6.55602482e-14,  4.33210907e-14,  2.14685336e-06],
           [-8.76525384e-02, -8.76525384e-02, -8.76533564e-02],
           [-8.76531009e-02,  8.76531009e-02,  8.76522830e-02],
           [ 8.76531009e-02, -8.76531009e-02,  8.76522830e-02],
           [ 8.76525384e-02,  8.76525384e-02, -8.76533564e-02]])

.. TODO: partial charges

Nice! We got our first results for `this methane molecule <https://cccbdb.nist.gov/exp2x.asp?casno=74828&charge=0>`_.
What more can we do?
Before we move to show the additional functionality we bring with :class:`Fireball`
we need to discuss some of the different options available.
In the next section we will start by explaining what *FData* are and the crucial role they play.
