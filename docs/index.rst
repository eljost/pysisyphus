.. pysisyphus documentation master file, created by
   sphinx-quickstart on Thu May 24 09:59:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pysisyphus documentation!
======================================

`pysisyphus` is a software-suite for the exploration of potential energy surfaces in ground-
and **excited states**. It implements several methods to search for stationary points
(minima and first order saddle points) and the calculation of minimum energy paths by means
of IRC and Chain of States methods like Nudged Elastic Band and Growing String.
Furthermore it provides tools to easily analyze & modify geometries (aligning, translating,
**interpolating**, ...) and to visualize the calculation results/progress.

The required energies, gradients and hessians are calculated by calling external quantum chemistry codes. `pysisyphus` can also be used as a library to implement custom quantum chemistry workflows.

If any issues arise please open an issue_ and I'll see if it can be fixed and my time permits it.
Contrubtions are welcome, so feel free to submit a PR.

**This software is still work in progress. Use at your own risk. Also take a look at the** `license <https://github.com/eljost/pysisyphus/blob/master/LICENSE>`_

.. _license: https://github.com/eljost/pysisyphus/blob/master/LICENSE
.. _issue: https://github.com/eljost/pysisyphus/issues

.. toctree::
   :maxdepth: 3
   :numbered:
   :caption: Contents:

   overview.rst
   installation.rst
   nix.rst
   calculators.rst
   min_optimization.rst
   es_optimization.rst
   chainofstates.rst
   tsoptimization.rst
   irc.rst
   plotting.rst
   pysistrj.rst

Indices and tables
==================

.. * :ref:`genindex`
* :ref:`search`
