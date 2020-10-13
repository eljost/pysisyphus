Chain Of States Methods
***********************

A chain of states (COS) describes a collection of distinct images
describing a chemical transformation.

When properly relaxed a COS may coincide with a minimum energy path (MEP), or is
(hopefully) a good approximation to it.
Tangents can be defined for every COS image and together they make up a discretized
path that describes the reaction/chemical transformation. The tangents are also used to
divide the COS gradient into components perpendicular and parallel to the discretized
path. Relaxation (optimization) of the COS is achieved by minimizing the perpendicular
component of the gradient, so only the parallel component remains.

With different Nudged Elastic Band (NEB) and the String Methods (SM) `pysisyphus` implements
several COS methods.

String method
=============

(When discussing String methods the COS 'images' are usually called 'nodes'.)

In the string method the whole COS is regularly reparametrized by fitting a spline
through all nodes and redistributing them along the spline according to a predefined
parametrization. By varying the parametrization equal node spacing or a higher
resolution around the highest energy image (HEI) can be achieved.
The tangents needed for the gradient projection are obtained as first derivatives
of the spline.

Reparametrization every :math:`n`-th cycle impedes the efficient optimization
of the string, this prevents the use of optimizers with some kind of history
like BFGS, L-BFGS or Conjugate Gradient (CG). The history has to be reset after
each reparametrization and only a simple Steepest Descent (SD) step can be done.

Nudeged Elastic Band
====================

No reparametrization takes place in the NEB method. It is characterized by the
replacement of the parallel gradient component by an articial spring force.
In principle optimizing a NEB should be easier as there is no reparametrization
and more sophisticated optimizers beyond SD can and should be employed.

General remarks
===============

Converged COS can yield good guesses for subsequent TS searches by using the HEI
as TS guess. In `pysisyphus` a subsequent TS search is easily started by including
`tsopt:` in the YAML input.
Knowledge of the initial and final images of the COS is used to construct a more
complete set of internal coordinates for the TS guess and it is less likely that
important coordinates are missed.
The initial imaginary mode to follow uphill is selected as the one featuring the
highest overlap with HEI tangent.

YAML example(s)
===============

Below you can find an example YAML-input including the most important options
that the user may want to modify when running a COS optimization.

.. code:: yaml

    preopt:                                  # Preoptimize inital and final geometry
    cos:
     type: gs                                # Do a growing string
     max_nodes: 9                            # Total string will have max_nodes + 2 images
    opt:
     type: string                            # Optimizer for GrowingString
     stop_in_when_full: 0                    # Stop string optimization N cycles after fully grown
     align: False                            # Disable Kabsch algorithm. Should be True with
                                             # coord_type == cartesian
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
    xyz: [first_preopt.xyz, last_preopt.xyz]
    coord_type: dlc                         # Run GrowingString in delocalized internal coordinates
                                            # (preferred).

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
