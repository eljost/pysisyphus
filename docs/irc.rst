Intrinsic Reaction Coordinate (IRC)
***********************************

By default two paths are integrated in plus- and minus-direction of the imaginary
mode (transition vector, TV) at the transition state (TS). If the IRC is started from
a non-stationary point with non-vanishing gradient the `downhill: True` argument can be
set to integrate only one path and skip the initial hessian calculation.

IRCs are integrated in mass-weighted coordinates and the default integrator is
`EulerPC`.

The initial displacement from the TS is done by requiring a certain energy lowering
:math:`\mathrm{d}E` from moving along the TV and calculating the corresponding step
length from a quadratic potential: :math:`\mathrm{d}E = \frac{1}{2} (k \cdot \mathrm{d}q^2)`,
with :math:`k` being the eigenvalue of the TV (imaginary mode) and :math:`\mathrm{d}q`
the required step length. The default required energy lowering is
:math:`0.0005 ~ \mathrm{au}`. Alternatively the initial can be done by a prescribed length:
`displ: length` and `displ_length: [desired step]` (default is :math:`0.1 \sqrt{\mathrm{amu}}
\cdot \mathrm{bohr}`).

The resulting endpoints of the IRC integration can be further optimized to stationary poins
by adding the `endopt:` section (vide infra). By setting `fragments: True` in `endopt` separate
fragments will be detected and optimized individually. This may be useful if the molecules are
only weakly bound.

By default IRC path data is dumped to `dump_fn: irc_data.h5` every `dump_every: 5` cycles.
IRC paths can be plotted with `pysisplot --irc`.

YAML example(s)
===============

Below you can find an example YAML-input including the most important options
that the user may want to modify.

.. code:: yaml

    irc:
     type: eulerpc                      # Similar to EulerPC from Gaussian

     #displ: energy|length              # How to do the initial displacement
     #displ_energy: 0.0005              # Energy lowering in au (Hartree)
     #displ_length: 0.1                 # Step length along the TV

     #forward: True
     #backward: True
     #downhill: False                   # Only integrate downhill

     #rms_grad_thresh: 0.0005           # Convergence threshold
    endopt:
     #fragments: False           # Optimize potential fragments separately
     coord_type: redund_v2       # Use different coordinates if desired
     trust_max: 0.1
     thresh: gau_tight
    calc:
     type: xtb                          # extended tight-binding calculator
     pal: 4
     charge: 0
     mult: 1
    geom:
     fn: hfabstraction_ts_opt_xtb.xyz   # Input coordinates

Further examples for IRC calculations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/irc>`_.

IRC base class
=====================

Base class for IRC integrators from which actual IRC integrators are derived.

.. automodule:: pysisyphus.irc.IRC
    :members:
    :undoc-members:

IRC Integrators
==========================

Damped-Velocity-Verlet integrator
----------------------------

.. automodule:: pysisyphus.irc.DampedVelocityVerlet
    :members:
    :undoc-members:
    :show-inheritance:

Euler integrator
----------------------------

Not recommended as it only produces reasonable results with very small
step sizes.

.. automodule:: pysisyphus.irc.Euler
    :members:
    :undoc-members:
    :show-inheritance:

Euler-Predictor-Corrector integrator
----------------------------

Recommended IRC integrator and default choice.

.. automodule:: pysisyphus.irc.EulerPC
    :members:
    :undoc-members:
    :show-inheritance:

Gonzales-Schlegel-2 integrator
----------------------------

.. automodule:: pysisyphus.irc.GonzalesSchlegel
    :members:
    :undoc-members:
    :show-inheritance:

Local-Quadratic-Approximation integrator
----------------------------

.. automodule:: pysisyphus.irc.LQA
    :members:
    :undoc-members:
    :show-inheritance:

Modified-Ishida-Morokuma-Komornicki integrator
----------------------------

Similar to the algorithm implemented by ORCA 4.

.. automodule:: pysisyphus.irc.IMKMod
    :members:
    :undoc-members:
    :show-inheritance:

Runge-Kutta-4 integrator
----------------------------

Not recommended, as it is very slow.

.. automodule:: pysisyphus.irc.RK4
    :members:
    :undoc-members:
    :show-inheritance:
