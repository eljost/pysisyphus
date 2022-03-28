Minimization
************

Searching for minimum energy geometries of molecules is preferably done using
second order methods that employ hessian information. Using the plain hessian
:math:`H` for step prediction through :math:`p=-H^{-1}g`, with :math:`g` being
the gradient vector, may yield erroneous uphill steps when :math:`H` has negative
eigenvalues.

This can be understood by calculating the step in the basis of the hessian eigenvectors
:math:`V`.

.. math::

    \begin{align}
        \tilde{H} &= V^T H V \\
        \tilde{g} &= V^T g \\
        \tilde{p_i} &= -\frac{\tilde{g}_i}{\tilde{H}_{ii}} \\
    \end{align}

:math:`\tilde{H}, \tilde{g}` and :math:`\tilde{p}` are transformed into the eigenvector
basis and the subscript :math:`i` indicates the component belongs to the :math:`i`-th
eigenvalue (-vector) of :math:`H`. As the gradient always points into the direction of
greater function values dividing it by a negative eigenvalues :math:`\tilde{H}_{ii}` will
lead to a step in uphill direction along the :math:`i`-th eigenvector.
The step in the original basis is obtained by a simple back-transformation:

.. math::

        p = V \tilde{p}

A step in downhill direction can be ensured by introducing an appropriate shift parameter
:math:`\lambda` that must be smaller than the smallest eigenvalue of :math:`H`:

.. math::

    \tilde{p_i} = -\frac{\tilde{g}_i}{\tilde{H}_{ii} - \lambda} \\
    

In `pysisyphus` the shift parameter :math:`\lambda` is obtained by the Rational Function
Optimization approach.

When the root-mean-square of the predicted step :math:`p` drops below a certain threshold
(default is 0.0025 au or rad) controlled GDIIS is tried. If GDIIS fails or is not yet
possible a polynomial line-search is conducted. Using the latest two energies and the
projection of the latest two gradients onto the last step the coefficients of a cubic
polynomial are defined unambiguously. By requiring :math:`\mathrm{d}^2 f(x)/\mathrm{d}x^2 >= 0`
and the equality holding at exactly one point also a quartic polynomial can be determined.

For now line-searches using a cubic polynomial are disabled as they seem to degrade the
optimizer performance, so only the constrained quartic polynomial is tried. By also using
projected hessian information the coefficients of a quintic polynomial can be determined,
but this also seems to be inferior compared to the constrained quartic polynomial.

Convergence is indicated when the root-mean-square and the maximum value of the gradient
and step drop below a certain threshold. It is not uncommon that the gradient convergence
is achieved before step convergence. By using `overachieve_factor: [n, float > 1]` in
the YAML input convergence will be signalled, when gradient convergence is overachieved
by factor `n`.

YAML Example
------------

Below you can find an example YAML-input including the most important options
that the user may want to modify when using the RFOptimizer.

.. code:: yaml

    opt:
     type: rfo                      # Optimization algorithm
     max_cycles: 50                 # Maximum number of optimization cycles

     overachieve_factor: 2          # Indicate convergence, regardless of the 
                                    # proposed step when max(grad) and rms(grad)
                                    # are overachieved by factor [n]

     do_hess: True                  # Calculate the hessian at the final geometry
                                    # after the optimization.

     #hessian_recalc: None          # Recalculate exact hessian every n-th cylce

     #hessian_recalc_adapt: None    # Expects a float. Recalculate exact hessian
                                    # whenever the gradient norms drops below
                                    # 1/[n] of the gradient norm at the last hessian
                                    # recalculation.

     #hessian_init: fischer         # Type of model hessian. Other options are: 'calc,
                                    # simple, xtb, lindh, swart, unit'

     #hessian_update: bfgs          # Hessian-update. Other options are: 'flowchart,
                                    # damped_bfgs, bofill'. bofill is not recommended
                                    # for minimum searches.
     #small_eigval_thresh: 1e-8     # Neglect eigenvalues and corresponding eigenvectors
                                    # below this threshold.

     #max_micro_cycles: 50          # No. of micro cycles for the RS-variants. Does not apply
                                    # to TRIM.

     #trust_radius: 0.3             # Initial trust radius.
     #trust_max: 1.0                # Max. trust radius
     #trust_min: 0.1                # Min. trust radius

     #line_search: True             # Do line search

     #gdiis_thresh: 0.0025          # May do GDIIS if rms(step) falls below this threshold
     #gediis_thresh: 0.01           # May do GEDIIS if rms(grad) falls below this threshold
     #gdiis: True                   # Do controlled GDIIS after 'gdiis_thresh' is reached
     #gediis: False                 # Do GEDIIS after 'gediis_thresh' is reached

    calc:
     type: turbomole
     control_path: control_path_pbe0_def2svp_s1     # Path to the prepared calculation
     track: True                                    # Activate excited state tracking
     ovlp_type: tden                                # Track with transition density matrix overlaps
     charge: 0
     mult: 1
     pal: 4
     mem: 2000

    geom:
     type: redund
     fn: cytosin.xyz

Further examples for optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/opt>`_.

Convergence Thresholds
----------------------
Optimization convergence is indicated when four critera are met. Each
cycle the root-mean-square and the absolute maximum of force and proposed
step vectors are checked. The default criteria in pysisyphus correspond
to the `opt=loose` in Gaussian. Different thresholds can easily be requested
by supplying the respective keyword in the YAML input.

Convergence is also indicated, when max(force) and rms(force) are overachieved,
by a certain factor (`overachieve_factor`, default=3), similar to ORCAs optimizer.
Further keywords, controlling the convergence check are found below.

.. code:: yaml

    opt:
     ...
     thresh: gau_loose           # See table below for different thresholds.
     overachieve_factor: 3
     rms_force: null             # Derive four critera from this value
     rms_force_only: False       # Only check rms(force)
     max_force_only: False       # Only check max(force)
     # Chain-of-States specific
     coord_diff_thresh: 0.01     # Terminate on insignificant coordinate change
     reparam_thresh: 0.001       # Terminate on insignificant coordinate change upon reparametrization
     ...

============ ========== =========== ========= ========= =======
Threshold    max(force) rms(forces) max(step) rms(step) Comment
============ ========== =========== ========= ========= =======
`gau_loose`  2.5e-3     1.7e-3      1.0e-2    6.7e-3    default
`gau`        4.5e-4     3.0e-4      1.8e-3    1.2e-3
`gau_tight`  1.5e-5     1.0e-5      6.0e-5    4.0e-5
`gau_vtight` 2.0e-6     1.0e-6      6.0e-6    4.0e-6
`baker`      3.0e-4     2.0e-4      3.0e-4    2.0e-4    energy difference and step are also checked
`never`      2.0e-6     1.0e-6      6.0e-6    4.0e-6    dummy thresholds; never converges
============ ========== =========== ========= ========= =======

Available Optimizers
---------------------------

Base class for optimizers that utilize Hessian information.

.. automodule:: pysisyphus.optimizers.Optimizer
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: pysisyphus.optimizers.HessianOptimizer
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: pysisyphus.optimizers.RFOptimizer
    :members:
    :undoc-members:
    :show-inheritance:
