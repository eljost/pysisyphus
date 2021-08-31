Chain Of States Methods
***********************

A chain of states (COS) comprises a set of distinct states (images) and is usually
spanned between two minima on a potential energy surface (PES).

When properly relaxed, a COS coincides with a minimum energy path (MEP), or is a good
approximation to it.
Tangents can be defined for every COS image and together they make up a discretized
path that describes the reaction/chemical transformation. The tangents are also used to
divide the COS gradient into perpendicular and parallel components, w.r.t the tangents.
Relaxation (optimization) of the COS is achieved by minimizing the perpendicular
component of the gradient, so only the parallel component remains.

`pysisyphus` offers different COS implementations, namely Nudged Elasic Band (NEB),
including its Adaptive and Free-End (and Free-End-Adaptive) modificiations and different
flavors of String methods (Growing String Method, GSM, and Simple Zero Temperature String, SZTS).
The GSM implementation is also available for with internal coordinates. 

String method
=============

(When discussing String methods the COS 'images' are usually called 'nodes'.)

In the string method the whole COS is periodically (every n-th cycle) reparametrized
by fitting a spline through all nodes and redistributing the nodes them along the spline,
according to a predefined parametrization.
By choosing between different parametrizations equal node spacing or higher
resolution around the highest energy image (HEI) can be achieved.
The tangents needed for the gradient projection are obtained as first derivatives
of the spline.

Reparametrization every :math:`n`-th cycle impedes efficient string optimization, and
prevents the use of optimizers with some kind of history Conjugate Gradient (CG).
The optimizer history is reset after each reparametrization and a simple Steepest Descent
(SD) step is done after reparametrization.

Nudeged Elastic Band
====================

No reparametrization takes place in the NEB method. The parallel gradient component
along the tangent is projected out and replaced by an artificial spring force.
In principle, optimizing NEBs should be easier as there is no reparametrization
and more sophisticated optimizers beyond SD can and should be employed.

General remarks
===============

Converged COS produce good guesses for subsequent TS searches, when the (splined)
HEI is determined. In `pysisyphus` subsequent TS searches are easily started by including
`tsopt:` in the YAML input.
Knowledge of the initial and final images of the COS is used to construct a more
complete set of internal coordinates for the TS guess and it is less likely that
important coordinates are missed.
The initial imaginary mode to follow uphill is selected as the one featuring the
highest overlap with HEI tangent.

YAML example(s)
===============

Below you can find an example YAML-input including the most important options
that the user may want to modify when running a GSM optimization.

.. code:: yaml

    precontr:                                # Preconditioning of translation & rotation
    preopt:                                  # Preoptimize inital and final geometry
    cos:
     type: gs                                # Do a growing string
     max_nodes: 9                            # Total string will have 9 + 2 == 11 images
     climb: False                            # Enable climbing image (CI), usually a good idea.
     climb_rms: 0.005                        # rms(forces) threshold for enabling CI
     climb_lanczos: False                    # Use tangent obtained from Lanczos algorithm for CI.
     climb_lanczos_rms: 0.005                # rms(forces) threshold for enabling Lanczos algorithm.
     reparam_check: rms                      # Criterian for growing new frontier nodes (rms/norm).
     perp_thresh: 0.05                       # Threshold for growing new frontier nodes.
     reparam_every: 2                        # Reparametrize every n-th cycle when NOT fully grown.
     reparam_every_full: 3                   # Reparametrize every n-th cycle when fully grown.
    opt:
     type: string                            # Optimizer for GrowingString
     stop_in_when_full: -1                   # Stop string optimization N cycles after fully grown
                                             # Disabled by -1. Usually it's a better idea to
                                             # further converge the string, after is fully grown.
     align: False                            # Disable Kabsch algorithm. Should be True with
                                             # type: cart
     scale_step: global                      # Scale down step as whole (global) or per image
                                             # (per_image)
    tsopt:
     type: rsprfo                            # Continue with TS-optimization of highest energy images
                                             # (HEI) using the RS-P-RFO algorithm
     do_hess: True                           # Calculate hessian at optimized TS geometry
     trust_max: 0.3
     thresh: gau_loose
    calc:
     type: orca
     keywords: "b3lyp 6-31G* rijcosx"
     pal: 4
     charge: 0
     mult: 1
    geom:
     type: dlc
     fn: [first_preopt.xyz, last_preopt.xyz] # Run GrowingString in delocalized internal coordinates
                                             # (preferred).

For NEB optimizations a different optimizer (not :code:`type: string`) should be used, e.g.,
QuickMin :code:`type: qm` in the beginning, and :code:`type: lbfgs` later on. It is not
yet possible to specify two different optimizers that are used in stages, so if is
desired it must be done manually.

.. code:: yaml

    # Taken from examples/complex/06_diels_alder...
    geom:
     type: cart
     fn: diels_alder.trj
    calc:
     type: xtb
     charge: 0
     mult: 1
     pal: 4
    preopt:
     max_cycles: 5
    interpol:                       # In NEBs the whole path is interpolated beforehand.
     type: redund                   # Possible values: redund|idpp|lst|linear
     between: 10                    # Interpolate n-geometries between every pair of supplied
                                    # geometries. For two this yields 1 + 10 + 1 == 12 images,
                                    # for three geometries this yields 1 + 10 + 1 + 10 + 1 == 23
                                    # geometries.
    cos:
     type: neb
     climb: False
    opt:
     type: lbfgs
     align: True
     rms_force: 0.01
     max_step: 0.04
    tsopt:
     type: rsirfo
     do_hess: True
     max_cycles: 75
     thresh: gau_tight
     hessian_recalc: 7
    irc:
     type: eulerpc
     rms_grad_thresh: 0.0005 
    endopt:

Further examples for COS optimizations from `.yaml` input can be found
`here <https://github.com/eljost/pysisyphus/tree/dev/examples/complex>`_.

General advice for COS optimizations
===================================

- Start from optimized geometries or use the `preopt:` key in the YAML input.
- Consider fixing the initial and final images `fix_ends: True`
- Always use `align: True` when optimizing a COS in cartesian coordinates to remove
  translation and rotation. `align: True` must not be used when running a COS with DLC.
- Don't over-converge a COS. It is usually a better idea to converge a COS loosely
  and use the highest energy image (HEI) as a guess for a subsequent transition state
  (TS) search.
- If possible use a climbing image (CI) `climb: True`
- When running a growing string calculation (`type: gs`) use `stop_in_when_full: [n]` in
  the `opt:` section with a small integer `[n]` to stop the COS relaxation after the string
  is fully grown

Chain Of States base class
==========================

Base class for chain of state methods

.. automodule:: pysisyphus.cos.ChainOfStates
    :members:
    :undoc-members:

Chain Of State Methods
======================

Nudged Elastic Band (NEB)
-------------------------

.. automodule:: pysisyphus.cos.NEB
    :members:
    :undoc-members:
    :show-inheritance:

Adaptive NEB
------------

.. automodule:: pysisyphus.cos.AdaptiveNEB
    :members:
    :undoc-members:
    :show-inheritance:

Free-End NEB
------------

.. automodule:: pysisyphus.cos.FreeEndNEB
    :members:
    :undoc-members:
    :show-inheritance:

Simple Zero-Temperature String
------------------------------

.. automodule:: pysisyphus.cos.SimpleZTS
    :members:
    :undoc-members:
    :show-inheritance:

Growing Chain Of States base class
==========================

Base class for growing chain of state methods

.. automodule:: pysisyphus.cos.GrowingChainOfStates
    :members:
    :undoc-members:

Growing Chain Of State Methods
======================

Growing String Method
-------------------------

.. automodule:: pysisyphus.cos.GrowingString
    :members:
    :undoc-members:
    :show-inheritance:
