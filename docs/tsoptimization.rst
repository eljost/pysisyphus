Transition State Optimization
*****************************

Optimizing transition states (TSs) is black magic, especially in internal
coordinates. The most promising TS optimizers employ hessian information,
but locating TSs in `pysisyphus` is also possible with only first derivatives.

YAML example(s)
==============

Below you can find an example YAML-input including the most important options
that the user may want to modify for the hessian-based optimizers (RS-I-RFO,
RS-P-RFO, TRIM).

.. code:: yaml

    tsopt:
     type: rsirfo|rsprfo|trim       # Optimization algorithm

     do_hess: True                  # Calculate the hessian at the final geometry
                                    # after the optimization.

     hessian_recalc: 1              # Recalculate the exact hessian every n-th cylce
                                    # Very costly but a small number like 5 may be a good
                                    # idea.

     prim_coord: [3, 18]            # Select the mode to follow uphill by overlap with
                                    # this primitive internal coordinate

     #root: 0                       # Follow the n-th imaginary mode uphill, minimize
                                    # along the others.

     #hessian_init: calc            # Initial hessian is calculated exactly.

     #hessian_update: bofill        # Bofill hessian-update. Should not be modified for
                                    # TS optimizations.

     #hessian_ref: [path]           # Expects a path to a precalculated hessian at any
                                    # level of theory. Mode to maximize along is selected
                                    # by highest overlap with imaginary mode from the reference
                                    # hessian.

     #rx_coords: [[3, 18],]         # Skip calculation of hessian and modify a model hessian
                                    # maximize along the given primitive internals. Has to be
                                    # used with 'hess_init: swart|fischer|lindh|xtb' to avoid
                                    # calculation of exact hessian. Also requires 'coord_type:
                                    # redund|dlc' to work.

     #max_micro_cycles: 50          # No. of micro cycles for the RS-variants. Does not apply
                                    # to TRIM.

     #trust_radius: 0.3             # Initial trust radius.
     #trust_max: 1.0                # Max. trust radius
     #trust_min: 0.1                # Min. trust radius
    calc:
     type: xtb
     pal: 4
     charge: 0
     mult: 1
    xyz: shaked.geom_000.xyz
    add_prims: [[24, 20], ]         # If using internal coordinates ALWAYS check the coordinates
                                    # that pysisyphus generates (pysistrj [xyz] --internals). If
                                    # some important (reaction) coordinates appears to be missing
                                    # define them manually.

    coord_type: redund              # Optimization in internal coordinates.

Further examples for TS optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/tsopt>`_.

TSHessianOptimizer base class
=============================

Base class for TS optimizers that empoly a hessian.

.. automodule:: pysisyphus.tsoptimizers.TSHessianOptimizer
    :members:
    :undoc-members:
    :show-inheritance:

TS-Optimizer using hessian information
======================================

Restricted-Step Partitioned-RFO
-------------------------------

.. automodule:: pysisyphus.tsoptimizers.RSPRFOptimizer
    :members:
    :undoc-members:
    :show-inheritance:

Restricted-Step Image-Function-RFO
----------------------------------

.. automodule:: pysisyphus.tsoptimizers.RSIRFOptimizer
    :members:
    :undoc-members:
    :show-inheritance:

Trust-Region Image-Function Method
----------------------------------

.. automodule:: pysisyphus.tsoptimizers.TRIM
    :members:
    :undoc-members:
    :show-inheritance:

TS-Optimizers without hessian information
=========================================

Dimer method
------------

.. automodule:: pysisyphus.tsoptimizers.dimer
    :members:
    :undoc-members:
