#!/usr/bin/env python3

from setuptools import find_packages, setup
import sys

import versioneer

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

setup(
    name="pysisyphus",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python implementation of NEB and IRC algorithms.",
    url="https://github.com/eljost/pysisyphus",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "attrs",
        "distributed",
        "h5py",
        "numpy",
        "matplotlib",
        "natsort",
        "dask",
        "pandas",
        "pytest",
        "pyyaml",
        "rmsd",
        "sympy",
        "scipy",
    ],
    entry_points={
        "console_scripts": [
            "pysis = pysisyphus.run:run",
            "pysisplot = pysisyphus.plot:run",
            "pysistrj = pysisyphus.trj:run",
            "pysisfilter = pysisyphus.filtertrj:run",
        ]
    },
)
