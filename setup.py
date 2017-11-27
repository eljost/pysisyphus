#!/usr/bin/env python3

from setuptools import find_packages, setup
import sys

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

setup(
    name="pysisyphus",
    version="0.0.5",
    description="Python implementation of NEB and IRC algorithms.",
    url="https://github.com/eljost/pysisyphus",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "scipy",
        "numpy",
        "pytest",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "pysisrun = pysisyphus.run:run",
            "pysisplot = pysisyphus.plot:run",
        ]
    },
)
