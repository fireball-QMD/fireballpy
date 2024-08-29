.. _fireballpy_docs_mixer:

****************
Mixer parameters
****************

Theoretical background
======================

In order to perform the self-consistent loop, the code needs to find the
point-charges per orbital which satisfy the fixed-point equation,

.. math::

   g(x) = x \Rightarrow F(x) \equiv g(x) - x = 0 \;,

where in our case :math:`x` and :math:`g(x)` will correspond to the
input and output shell charges vector, respectively.
This problem may be seen as a root-finding problem for :math:`F(x)`,
the difference between output and input shell charges.

There are many algorithms for solving this kind of problems, with Broyden
methods [Broyden1965]_ as the most well-known.
Nonetheless, for charge mixing it is found that those simple iterative methods
are not enough as keeping information from previous iterations accelerates and
provides additional stability [Vanderbilt1984]_.

The implemented algorithm in FireballPy for charge mixing is described in detail
by Johnson [Johnson1988]_.
As a brief overview, we implement the following equations,

.. math::

   \lvert Q^{(k+1)} \rangle &= \lvert Q^{(k)} \rangle + G^{(1)} \lvert F^{(k)} \rangle - \sum_{\max(1,k-m+1)\leq i\leq k-1} w_i \langle \alpha^{(k)}_i \,\vert\, u^{(i)} \rangle \,, \\
   \lvert u^{(i)} \rangle &= G^{(1)} \lvert \Delta F^{(i)} \rangle + \lvert \Delta Q^{(i)} \rangle \,, \\
   \lvert \Delta Q^{(i)} \rangle &= \lvert Q^{(i)} \rangle - \lvert Q^{(i-1)} \rangle \,, \\
   \lvert \Delta F^{(i)} \rangle &= \lvert F^{(i)} \rangle - \lvert F^{(i-1)} \rangle \,,


where :math:`m` is the mixing order, :math:`Q^{(k)}` is the :math:`k^{th}` iteration input shell charges vector,
:math:`F^{(k)} = F\left(Q^{(k)}\right)` is the difference between output and input shell charges vector at the :math:`k^{th}` iteration
and :math:`\alpha^{(k)}_i` the :math:`i^{th}` element of the solution of the following linear system of equations:

.. math::

   A\mathbf{f}^{(k)} &= \mathbf{\alpha}^{(k)} \,, \\
   A_{ij} &= w_0^2\delta_{ij} + w_i w_j \langle \Delta F^{(i)} \, \vert \, \Delta F^{(j)} \rangle \,, \\
   f^{(k)}_i &= w_i \langle \Delta F^{(i)} \, \vert \, F^{(k)} \rangle \,,

which needs to be solved at each iteration.

The final step is choosing appropriate initial guess of the inverse Jacobian matrix, :math:`G^{(1)}`,
the weights :math:`w_0` and :math:`w_i`, and the mixing order, :math:`m`.

For efficiency reasons it is interesting to choose :math:`G^{(1)} = \beta I`.
This :math:`\beta` controls how fast the algorithm will converge to the solution, but high values may lead to instabilities.
We have found after exhaustive testing that a value of 0.1 results in a good balance between stability on challenging systems and speed.

For the weights we found that :math:`w_0 = 0.01`, :math:`w_i = \langle \Delta F^{(i)} \, \vert \, \Delta F^{(i)} \rangle^{-1/2}`,
provide with very fast convergence.
Nonetheless, this is not the only available choice.
For example, fixing :math:`w_0 = 0` and :math:`w_i = 1` will result in the well-known Anderson acceleration.
Also, :math:`w_0 = 1` and :math:`w_i = 0` leads to a suboptimal version of the Louie algorithm, which may suffice for very broad approximations.

Lastly, the for the mixing order we declare that 6 is optimal in both speed and stability.
However, some authors prefer 3 or even 2, which may be a good choice for highly difficult systems.

.. [Broyden1965] C, G, Broyden (1965). *A class of methods for solving
                 nonlinear simultaneous equations.*
                 `doi:10.1090/S0025-5718-1965-0198670-6 <https://doi.org/10.1090/S0025-5718-1965-0198670-6>_`

.. [Vanderbilt1984] David Vanderbilt and Steven G. Louie (1984).
                    *Total energies of diamond (111) surface
                    reconstructions by a linear combination of
                    atomic orbital method.*
                    `doi:10.1103/PhysRevB.30.6118 <https://doi.org/10.1103/PhysRevB.30.6118>`_

.. [Johnson1988] D.D. Johnson (1988). *Modified Broyden's method
                 for accelerating convergence in self-consistent
                 calculations.*
                 `doi:10.1103/PhysRevB.38.12807 <https://doi.org/10.1103/PhysRevB.38.12807>`_

Tweaking Mixer Parameters
=========================

Mixer parameters in FireballPy are controlled by the
dictionary :code:`mixer_kws`.
It has the following keys:

- :code:`method`: :code:`string` Can be either 'anderson' or 'johnson'.
  This will automatically set the weights to replicate
  Anderson acceleration ('anderson') or
  ours based on Johnson's findings ('johnson').
  If the :code:`mixer_kws` input is not provided, it will assume Johnson standards.

- :code:`mix_order`: :code:`int` Number of previous steps to mix
  in the solution (:math:`m`). Default is 6.

- :code:`beta`: :code:`float` Scale factor of the identity matrix
  as initial guess for the inverse Jacobian matrix (:math:`\beta`).
  Default is 0.1.

- :code:`tol`: :code:`float` Convergence tolerance to stop the algorithm.
  Default is :math:`10^{-8}`.

- :code:`max_iter`: :code:`int` Maximum number of iterations allowed if convergence is not reached.
  Default is 200.

- :code:`w0`: :code:`float` Value of the zeroth weight (:math:`w_0`).
  This option will be ignored if :code:`method` is set.





























