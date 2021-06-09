# pysisyphus
[![Documentation Status](https://readthedocs.org/projects/pysisyphus/badge/?version=master)](https://pysisyphus.readthedocs.io/en/master/?badge=master)
[![build](https://github.com/eljost/pysisyphus/workflows/Python%20application/badge.svg)](https://github.com/eljost/pysisyphus/actions)
[![codecov](https://codecov.io/gh/eljost/pysisyphus/branch/master/graph/badge.svg)](https://codecov.io/gh/eljost/pysisyphus)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/96281078.svg)](https://zenodo.org/badge/latestdoi/96281078)
[![Join the chat at https://gitter.im/eljost/pysisyphus](https://badges.gitter.im/eljost/pysisyphus.svg)](https://gitter.im/eljost/pysisyphus?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

`pysisyphus` is a software-suite for the exploration of potential energy surfaces in ground-
and **excited states**. It implements several methods to search for stationary points
(minima and first order saddle points) and the calculation of minimum energy paths by means
of IRC and Chain of States methods like Nudged Elastic Band and Growing String.
Furthermore it provides tools to easily analyze & modify geometries (aligning, translating, **interpolating**, ...) and to visualize the calculation results/progress.

**Further information can be found in the Open Access [pysisyphus paper](https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26390).**

Required energies, gradients and hessians are calculated by calling external quantum chemistry codes. Alternatively `pysisyphus` can also be used as a library to implement custom quantum chemistry workflows.

If any issues arise please open an issue and I'll try to fix it if possible and my time permits it. Contrubtions are welcome, so feel free to submit a PR.

**This software is still work in progress. Use at your own risk.**

## Example

Fully internal coordinate transition state search for the famous alanine dipeptide isomerization reaction, using the xtb calculator and the growing string method.

[![asciicast](https://asciinema.org/a/300731.svg)](https://asciinema.org/a/300731)

## Documentation

Please have a look at the [documentation](https://pysisyphus.readthedocs.io/en/dev/).
