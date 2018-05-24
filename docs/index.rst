.. pysisyphus documentation master file, created by
   sphinx-quickstart on Thu May 24 09:59:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pysisyphus's documentation!
======================================

pysisphus implements Chain Of States (COS) methods like Nudged Elastic Band
(NEB) and Simple Zero Temperature String (SZTS) to converge minimum energy
paths and provide initial guesses for transition states (TS) optimizations.

In addition pysisyphus provides serveral Intrinsic Reaction Coordinate
algorithms. The required gradients and/or hessians are calculated by
calling external quantum chemistry codes. By default everything is done
in cartesian coordinates but an internal coordinates implementation is
in progress.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
 
   installation.rst
   calculators.rst
   optimizers.rst
   examples.rst
   pysisyphus
	


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
