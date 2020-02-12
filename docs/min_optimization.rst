Optimization of Minima
**********************

Searching for minimum energy geometries of molecules is preferably done using
second order methods that employ hessian information. Using the plain hessian
:math:`H` for step prediction through :math:`p=-H^{-1}g`, with :math:`g` being
the gradient, may yield erroneous uphill steps when :math:`H` has negative eigenvalues.

This can be understood by calculating the step in the basis of the hessian eigenvectors
:math:`V`.

.. math::

    \begin{align}
        \tilde{H} &= V^T H V \\
        \tilde{g} &= V^T g \\
        \tilde{p_i} &= -\frac{\tilde{g}_i}{\tilde{H}_{ii}} \\
    \end{align}

:math:`\tilde{H}, \tilde{g}` and :math:`\tilde{p}` are transformed into the eigenvector
basis and the subscript :math:`_i` indicates the component belongs to the :math:`i`-th
eigenvalue (-vector) of :math:`H`. As the gradient always points into the direction of
greater function values dividing it by negative eigenvalues :math:`\tilde{H}_{ii}` will
lead to a step in uphill direction along the :math:`i`-th eigenvector of :math:`H`.

The step in the original basis is obtained by a simple back-transformation:

.. math::

        p = V \tilde{p}

A step in downhill direction can be ensured by introducing a shift parameter :math:`\lambda`
that must smaller than the smallest eigenvalue :math:`H_{ii}`:

.. math::

    \tilde{p_i} = -\frac{\tilde{g}_i}{\tilde{H}_{ii} - \lambda} \\
    

In `pysisyphus` shift parameter :math:`\lambda` is obtained from the Rational Function
Optimization approach.


YAML example
===============

Below you can find an example YAML-input including the most important options
that the user may want to modify when using the RFOptimizer.

.. code:: yaml

    opt:
     type: rfo                      # Optimization algorithm

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

    xyz: cytosin.xyz
    coord_type: redund

Further examples for optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/opt>`_.

HessianOptimizer base class
=============================

Base class for optimizers that empoly a hessian.

.. automodule:: pysisyphus.optimizers.HessianOptimizer
    :members:
    :undoc-members:
    :show-inheritance:

Optimizer using hessian information
======================================

Restricted-Step RFO
-------------------------------

.. automodule:: pysisyphus.optimizers.RFOptimizer
    :members:
    :undoc-members:
    :show-inheritance:
