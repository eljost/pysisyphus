#!/usr/bin/env python3

from os import path
from setuptools import find_packages, setup
import sys

import versioneer

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

# Read contents of the README for use in the description on PyPI
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as handle:
    long_description = handle.read()

setup(
    name="pysisyphus",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python suite for exploring potential energy surfaces.",
    long_description=long_description,
    long_description_content_type='text/markdown',
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
