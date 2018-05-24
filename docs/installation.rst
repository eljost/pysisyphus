Installation
===========

Dependencies
------------

::

	scipy
	numpy
	pytest
	matplotlib
	natsort
	pyyaml

Getting the code
-------------------

The **recommended way**, as the version on pypi may be outdated:

::

	cd [root dir where you want to keep pysisyphus]
	git clone https://github.com/eljost/pysisyphus
	cd pysisyphus
	python setup.py develop

or **(not recommended)**:

::

	pip install pysisyphus

Post-install configuration
--------------------------------------------

Create a file called ``.pysisyphusrc`` in your $HOME-directory
containing the commands to call your installed calculators, e.g.:

::

	[orca]
	cmd=/scratch/programme/orca-4.0.1.2/orca

	[xtb]
	cmd=xtb

	[openmolcas]
	cmd=pymolcas

	[gaussian09]
	cmd=/usr/local/gaussian/g09/g09

	[gaussian16]
	cmd=/usr/remote/gaussian/g16/g16

	[wfoverlap]
	cmd=/scratch/wfoverlap_1.0/bin/wfoverlap.x

When the specified path is not absolute but relative (as above for xtb
and openmolcas) the corresponding binaries have to be
available on the $PATH! The **Turbomole** binaries (dscf, ridft, ...) are
hardcoded in the pysisyphus calculator and don't have to be specified here,
but they still have to be available on $PATH.
