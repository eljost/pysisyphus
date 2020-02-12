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
    license="License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "attrs",
        "autograd",
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
    # Install locally with
    #   pip install -e .[extra]
    extras_require={
        "qcengine": ["qcengine>=0.13.0", ],
        "ase": ["ase", ],
        "pyscf": ["pyscf>=1.7.0", ],
    },
    entry_points={
        "console_scripts": [
            "pysis = pysisyphus.run:run",
            "pysisplot = pysisyphus.plot:run",
            "pysistrj = pysisyphus.trj:run",
            "pysisfilter = pysisyphus.filtertrj:run",
            "pysisserver = pysisyphus.server.main:run",
        ]
    },
)
