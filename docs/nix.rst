Installation via Nix
********************

An easy and reliable way to get a complete pysisyphus installation, including additional quantum chemistry (QC) codes, is possible through the `Nix package manager`_. Required Nix files for the QC codes and their dependencies (linear algebra & MPI libraries etc.) are provided by the NixOS-QChem_ overlay.

Nix installations are fully reproducible and extremely simple to accomplish for any potential pysisyphus user, as most of the necessary configuration is already described in the respective Nix files.

The `Nix Pills`_ series provides an informative introduction to the general concepts behind Nix.

Recent stable versions of pysisyphus are also directly included in NixOS-QChem_.
If you do not require a development version of pysisyphus, using the NixOS-QChem_ overlay is the preferred option.
Otherwise, continue with the instructions below.


Prerequisites
=============

Please follow the installation instructions in the `nix manual`_ to set up a working Nix installation on your system.
You may want to enable the :code:`flakes` and :code:`nix-command` features for Nix.

Caching
=======

The build time can be reduced drastically, if the quantum chemistry codes are fetched from a binary cache.

Pysisyphus provides a binary cache on Cachix.

Nix >= 2.6.0 with flakes directly supports the binary cache without further action required.

Alternatively, if you are allowed to add binary caches in nix, you may simply execute:

.. code-block:: bash

    nix run nixpkgs#cachix use pysisyphus

Usage
=====

You may use pysisyphus directly via Flakes without cloning any repository or local installation.
Some examples how to use pysisyphus with nix:

Inspect the structure of pysisyphus' flake:

.. code-block:: bash

    nix flake show github:eljost/pysisyphus

Directly run the :code:`pysis` command on an input :

.. code-block:: bash

  nix run github:eljost/pysisyphus input.yml

Start an interactive shell with a pysisyphus installation

.. code-block:: bash

  nix shell github:eljost/pysisyphus

For developing and hacking pysisyphus and for legacy Nix commands, first lone the pysisyphus repository

.. code-block:: bash

    git clone https://github.com/eljost/pysisyphus.git && cd pysisyphus

You may obtain a development shell for pysisyphus via

.. code-block:: bash

    nix develop

build pysisyphus in the legacy nix style

.. code-block:: bash

    nix-build ./nix/default.nix

or obtain a development legacy-style dev-shell

.. code-block:: bash

    nix-shell ./nix/shell.nix

The Flake also offers options to build pysisyphus with proprietary components such as ORCA or Turbomole, e.g.

.. code-block:: bash

    nix run .#pysisyphusOrca

or build singularity or docker containers:

.. code-block:: bash

    nix build .#pysisyphusSingularity


.. _`Nix package manager`: https://nixos.org/download.html
.. _`NixOS-QChem`: https://github.com/markuskowa/NixOS-QChem
.. _`nix-shell`: https://nixos.org/nix/manual/#sec-nix-shell
.. _`nix manual`: https://nixos.org/manual/nix/stable/
.. _`Nix Pills`: https://nixos.org/guides/nix-pills/index.html
.. _`Nix Bundle`: https://github.com/matthewbauer/nix-bundle
