# pysisphus
[![Documentation Status](https://readthedocs.org/projects/pysisyphus/badge/?version=master)](https://pysisyphus.readthedocs.io/en/master/?badge=master)
![build](https://github.com/eljost/pysisyphus/workflows/Python%20application/badge.svg)
[![codecov](https://codecov.io/gh/eljost/pysisyphus/branch/master/graph/badge.svg)](https://codecov.io/gh/eljost/pysisyphus)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/96281078.svg)](https://zenodo.org/badge/latestdoi/96281078)

`pysisyphus` is a software-suite for the exploration of potential energy surfaces in ground-
and **excited states**. It implements several methods to search for stationary points
(minima and first order saddle points) and calculation of minimum energy paths by means
of IRC and Chain of States methods like Nudged Elastic Band and Growing String.
Furthermore it provides tools to easily analyze & modify geometries (aligning, translating, interpolating, ...) and to visualize the calculation results/progress.

The required energies, gradients and hessians are calculated by calling external quantum chemistry codes. `pysisyphus` can also be used as a library to implement custom quantum chemistry workflows.

If any issues arise please open an issue_ and I'll try to fix it if possible and my time permits it. Contrubtions are welcome, so feel free to submit a PR.

**This software is still work in progress. Use at your own risk. Also take a look at the [license](https://github.com/eljost/pysisyphus/blob/master/LICENSE)**

Please have a look at the [documentation](https://pysisyphus.readthedocs.io/en/dev/).
