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
    include_package_data=True,
    install_requires=[
        "autograd",
        "dask",
        "distributed",
        "h5py==3.2.1",
        "jinja2",
        "matplotlib",
        "numpy>=1.18.1",
        "natsort",
        "pytest",
        "pyyaml",
        "rmsd",
        "scipy>=1.4.1",
        "sympy>=1.5.1",
        "scikit-learn>=0.24.2",
    ],
    # Install locally with
    #   pip install -e .[extra]
    extras_require={
        "qcengine": ["qcengine>=0.17.0", ],
        "ase": ["ase>=3.21.0", ],
        "pyscf": ["pyscf>=1.7.6", ],
        # If you want to build the documentation
        "sphinx": ["sphinx", "sphinx-rtd-theme"],
    },
    entry_points={
        "console_scripts": [
            "pysis = pysisyphus.run:run",
            "pysisplot = pysisyphus.plot:run",
            "pysistrj = pysisyphus.trj:run",
            "pysisfilter = pysisyphus.filtertrj:run",
            "pysispack = pysisyphus.pack:run",
            "pysisthermo = pysisyphus.drivers.thermo:run_thermo",
        ]
    },
)
