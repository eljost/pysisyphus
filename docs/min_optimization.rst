Optimization of Minima
**********************

To be done.

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
