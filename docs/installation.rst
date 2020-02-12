Installation
************

Preparing an environment
========================

It is good idea to install `pysisyphus` into a separate python environment,
whether it is an `Anaconda <https://www.anaconda.com/>`_ environment or a
`virtualenv <https://docs.python.org/3/library/venv.html>`_ environment.

.. code-block:: bash

    # Optional: Create separate Anaconda environment
    conda create -n pysis-env python=3.7
    activate pysis-env
    # or virtual environment
    python3 -m venv pysis-env
    source pysis-env/bin/activate

Installation from PyPI with pip
===============================

This installs the latest stable release as published on PyPI. If you don't want to
do any development this is probably the preferred way of installation.

.. code-block:: bash

    pip install pysisyphus
    # Installation of extras is also possible
    # pip install pysisyphus[qcengine]

Installation from source
========================

Decide on an **$install_dir** and clone the repository from github. If you want to change
the code after installation do an editable (-e) installation with pip.

.. code-block:: bash

    git clone https://github.com/eljost/pysisyphus.git $install_dir
    cd $install_dir
    # Install with -e if you want an editable installation
    pip install [-e] .
    # Installation of extras is also possible
    # pip install [-e] .[qcengine]

Setting up .pysisyphusrc
========================

`pysisyphus` interfaces several quantum chemistry codes and related software (Multiwfn, Jmol).
Software available to `pysisyphus` is registered in **$HOME/.pysisyphusrc**, so `pysisyphus`
knows which command to execute an/or where to find the binaries.

Depending on the software different choices were made how it is registered. An example **.pysisyphusrc** is given below, with a short comment for each software.

.. code-block:: ini

    # Excited state calculators

    [gaussian09]
    # Cmd to execute. Please ensure that g09 is on your $PATH.
    cmd=g09
    formchk_cmd=formchk
    unfchk_cmd=unfchk

    [gaussian16]
    # Cmds to execute. Please ensure that the binaries are found in your $PATH.
    cmd=g16
    formchk_cmd=formchk
    unfchk_cmd=unfchk

    [openmolcas]
    # Cmd to execute. Please ensure that pymolcas is on your $PATH.
    cmd=pymolcas

    [orca]
    # ORCA needs the full path to its binary, so please provide the full path.
    cmd=/scratch/programme/orca_4_2_0_linux_x86-64_openmpi314/orca

    #[pyscf]
    # pyscf must not have an explicit entry in the .pysisyphusrc. pysisyphus uses
    # the python API of pyscf, so it is mandatory that is installed in the same environment
    # as pysisyphus.

    #[turbomole]
    # Turbomole must not have an explicit entry in the .pysisyphusrc. The user has to take
    # care that everything is set up correctly, e.g. TM-binaries are on the PATH etc...
    # The respective commands are hardcoded into pysisyphus (dscf, ridft, ricc2, ...)

    # Ground state calculators

    [mopac]
    # Similar to Psi4. An example is given below.
    cmd=/user/johannes/bin/runmopac.sh

    [psi4]
    # As the Psi4 installation without conda is, to put it slightly, tricky it was
    # decided to allow the installation of Psi4 into a separate conda environment.
    # pysisyphus then creates a Psi4 input and sends it to the (bash)-script given below
    # that accepts/expects one argument. It is the responsibility of the scrip to activate
    # the appropriate conda environment and submit the Psi4 input. An example runpsi4.sh
    # script is given below.
    cmd=/user/johannes/bin/runpsi4.sh

    #[qcengine]
    # QCEngine must not have an entry explicit entry in the .pysisyphusrc. It is used
    # via its python interface and can be installed as an extra with pip (see above).
    # The user is reference to the QCEngine for any further questions.

    [xtb]
    # Cmd to execute. Please ensure that xtb is on your $PATH.
    cmd=xtb

    # Utilities

    [wfoverlap]
    # Cmd to execute. Please ensure that wfoverlap is on your $PATH. The binary/source
    # can be obtained from https://github.com/sharc-md/sharc/tree/master/bin
    cmd=/scratch/wfoverlap_1.0/bin/wfoverlap.x

    [multiwfn]
    # Cmd to execute. Please ensure that Multiwfn is on your $PATH.
    cmd=Multiwfn

    [jmol]
    # Cmd to execute. Please ensure that jmol is on your $PATH.
    cmd=jmol


When the specified path/cmd is not absolute but relative (e.g. for xtb, g16, ...) the corresponding
binaries have to be available on the **$PATH** and all other environment variables have to
be set up correctly by the user.

Example runpsi4.sh
==================

.. code-block:: bash

    #!/bin/bash

    # Afaik this doesn't work in non-interactive shells ...
    # See https://github.com/conda/conda/issues/8072
    # conda activate psi4
    source /scratch/programme/anaconda3/bin/activate psi4
    #conda activate psi4
    psi4 -o stdout $1

Example runmopac.sh
==================

.. code-block:: bash

    #!/bin/bash

    module purge
    module load mopac

    MOPAC2016.exe $1
