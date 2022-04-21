Transition State Optimization
*****************************

To cite Frank Jensen: "Locating transition states
(TSs) is black magic, especially in internal coordinates" and I can say this is true.
The most promising TS optimizers in `pysisyphus` employ second derivative information
(hessian) but locating TS is also possible using only first derivatives by means of the
`dimer method` (DM).

TS optimizations should preferably be carried out in internal coordinates (`type: redund|dlc`).
Before starting the actual TS search, one should **always** check the internal coordinates
that pysisyphus sets up, by executing `pysistrj [xyz] --internals`, with `[xyz]` corresponding
to the geometry you plan to use for the optimization. **Check carefully** if all supposed
reaction coordinates are present. **If they are missing** define them manually with the
`add_prims` key in the YAML input.

Calculation of the exact hessian can be avoided by using the DM by using
`rx_coords: [list of primitives]` in combination with `hessian_init: fischer|lindh|swart|simple`.
By using both options, a diagonal model hessian is calcualted and modified for use in
a TS optimization.

The DM can be used with a powerful preconditioned LBFGS optimizer (:code:`type: lbfgs`).

Hessian Based TS Optimization
=============================

YAML example
------------

Below you can find an example YAML-input including the most important options
that the user may want to modify for the Hessian-based optimizers (RS-I-RFO,
RS-P-RFO, TRIM).

.. code:: yaml

    tsopt:
     type: rsirfo|rsprfo|trim       # Optimization algorithm

     do_hess: True                  # Calculate the hessian at the final geometry
                                    # after the optimization.

     hessian_recalc: 1              # Recalculate the exact hessian every n-th cylce
                                    # Very costly but a small number like 5 may be a good
                                    # idea.

     prim_coord: [BOND, 3, 18]      # Select the mode to follow uphill by overlap with
                                    # this primitive internal coordinate. Expects exactly
                                    # one typed primitive.

    rx_mode: [[[BOND, 3, 18], 1]]   # Select initial mode based on overlaps with a
                                    # constructed mode. Can be seen as generalization of
                                    # prim_coord. Expects a list of pairs, with each pair
                                    # comprising a typed primitive and a phase factor.
                                    #
                                    # See examples/{05,06}... for examples.

     #rx_coords: [[BOND, 3, 18],]   # Obtain initial mode by modifying a model Hessian.  Has
                                    # to be used with 'hess_init: swart|fischer|lindh|xtb'
                                    # to avoid calculation of exact hessian.
                                    # Also requires 'redund' or 'dlc' to work.

     #root: 0                       # Follow the n-th imaginary mode uphill, minimize
                                    # along the others.

     #hessian_init: calc            # Initial hessian is calculated exactly.

     #hessian_update: bofill        # Bofill hessian-update. Should not be modified for
                                    # TS optimizations.

     #hessian_ref: [path]           # Expects a path to a precalculated hessian at any
                                    # level of theory. The emode to maximize along is selected
                                    # by highest overlap with imaginary mode from the reference
                                    # hessian.


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
    geom:
     type: redund
     fn: shaked.geom_000.xyz
     add_prims: [[24, 20], ]        # If using internal coordinates ALWAYS check the coordinates
                                    # that pysisyphus generates (pysistrj [xyz] --internals). If
                                    # some important (reaction) coordinates appears to be missing
                                    # define them manually. In this example
                                    # we add a bond. If additional internal coordinates can
                                    # derived from the added primitives pysisyphus will do
                                    # it.

Second-derivative-free TS optimization can be done using the `Dimer` method. Attached
you can find a sample input.


Further examples for TS optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/tsopt>`_.


Dimer Method
============
`pysisyphus` implements the dimer method (DM) as calculator, wrapping another calculator
for the actual energy and gradient calculations. All DM configurations have to be done
in the :code:`calc:` section. The best results are obtained by optimizing the DM with
LBFGS in connection with a preconditioner (PLBFGS), estimated from the Lindh model Hessian.
in constrast to the TS optimizers mentioned above, PBLFGS is configured in the :code:`opt`
section of the YAML input. DM related information is logged to `dimer.log`. The current
DM orientation is saved in files with extension `.N`. Animated `.trj` files of the DM
orientation are saved to `.N.trj`.

YAML example
------------
Below you can find an example for the DM, applied to the isomerization of HCN.
Default values are commented out. See `examples/tsopt/02_hcn_tsopt_dimer <https://github.com/eljost/pysisyphus/tree/master/examples/tsopt/02_hcn_tsopt_dimer>`_ for the full example.


.. code:: yaml

    opt:
     type: plbfgs    # Preconditioned LBFGS
     thresh: baker   # Baker threshold, comparable to the Gaussian defaults
     do_hess: True   # Calculate Hessian at the end. May not be applicable to systems
                     # where Hessian calculation is infeasible.
    calc:
     type: dimer     # Dimer calculator wrapping PySCF at HF/3-21G level of theory
     calc:
      type: pyscf
      basis: 321g
      pal: 2
      charge: 0
      mult: 1
     #N_raw: [file]                # Path to a file containing an initial dimer orientation
     #length: 0.0189               # Separation of dimer images from common midpoint
     #rotation_max_cycles: 15      # Maximum number of rotations
     #rotation_method: fourier
     #rotation_thresh: 1e-4        # rms(rotationa force) threshold
     #rotation_tol: 1              # Rotation tolerance in degree. If a proposed trial
                                   # rotation angle is below the tolerance the rotation
                                   # will be skipped.
     #rotation_max_element: 0.001  # Max. step element for "rotation_method: direct"
     #rotation_interpolate: True   # Interpolate force on image after rotation
     #rotation_remove_trans: True  # Remove overall translation from N-vector
     #seed: 20182503                # Seed for the RNG for reproducability reasons.
    geom:
     type: cart
     fn: 01_hcn.xyz

General advice for TS optimizations
===================================

- Use as many exact hessians as your computational budget allows (`hessian_recalc: [n]`)
- Use tight convergence settings for your SCF. In contrast to codes like ORCA `pysisyphus`
  does not enforces this automatically
- Grid-based methods (DFT, some RI-types) may need finer grids to reduce numerical noise
- Try to preoptimize the TS using a cheap method (xtb) that allows `hessian_recalc: 1`
- Decrease current & allowed trust radius (`trust_radius: 0.1` and `trust_max: 0.3|0.5`)
- Check for orrect coordinate setup

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

Dimer method
----------------

.. automodule:: pysisyphus.calculators.Dimer
    :members:
    :undoc-members:
