Transition State Optimization
*****************************

To cite Frank Jensen: "Locating transition states
(TSs) is black magic, especially in internal coordinates" and I can say this is true.
The most promising TS optimizers in `pysisyphus` employ second derivative information
(hessian) but locating TS is also possible using only first derivative information
by means of the `Dimer` method.

TS optimizations should preferably be done in internal coordinates (`coord_type: redund`).
Before starting the actual TS search one should **always** check the internal coordinates
that pysisyphus set up by running `pysistrj [xyz] --internals`, with `[xyz]` corresponding
to the TS guess for the optimization. **Check carefully** if all supposed reaction coordinates
are present. **If they are missing** define them manually with the `add_prims` key
in the YAML input.

Calculation of the exact hessian can be avoided by using `type: dimer` or by using
`rx_coords: [list of primitives]` in combination with `hessian_init: fischer|lindh|swart|simple`.
With the latter options a diagonal model hessian is set up and the signs of the entries
corresponding to the reaction coordinates are inverted.

Right now the `Dimer` is only tested with `coord_type: cart`.

YAML example(s)
===============

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
                                    # this primitive internal coordinate. Expects indices
                                    # for one primitive internal.

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
                                    # define them manually. In this example
                                    # we add a bond. If additional internal coordinates can
                                    # derived from the added primitives pysisyphus will do
                                    # it.

    coord_type: redund              # Optimization in internal coordinates.

Second-derivative-free TS optimization can be done using the `Dimer` method. Attached
you can find a sample input.


Further examples for TS optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/tsopt>`_.

General advice for TS optimizations
===================================

- Use as many exact hessians as your computational budget allows it (`hessian_recalc: [n]`)
- Use tight convergence settings for your SCF. In contrast to codes like ORCA `pysisyphus`
  does not do this automatically
- Grid-based methods (DFT, some RI-types) may need finer grids to reduce numerical noice
- Try to preoptimize the TS using a cheap method (xtb) that allows `hessian_recalc: 1`
- Maybe decrease current & allowed trust radius (`trust_radius: 0.1` and `trust_max: 0.3|0.5`)
- Check for a correct coordinate setup

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

TS-Optimization without hessian information
===========================================

Old Dimer method
----------------

.. automodule:: pysisyphus.tsoptimizers.dimer
    :members:
    :undoc-members:

New Dimer method
----------------

.. automodule:: pysisyphus.calculators.Dimer
    :members:
    :undoc-members:
