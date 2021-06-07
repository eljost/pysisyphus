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

Configuration
=============

Typically no, or only very few, adjustements are necessary.

- :code:`optAVX` can be used to enable/disable AVX tunings in underlying quantum chemistry packages.
- :code:`licMolpro` can be set to a molpro license token to enable the Molpro package.
- :code:`optpath` can be set and point to the top-level installation directory of a Gaussian installation.

For more details on the configuration options see the NixOS-QChem_ repository.

The :code:`postOverlays` attribute can be used to further customise the underlying package sets.


Installation
============

The initial build can take some time, as basic computational chemistry packages will be downloaded and build.

Clone the pysisyphus repository

.. code-block:: bash

    git clone https://github.com/eljost/pysisyphus.git

and build the nix-derivation

.. code-block:: bash

    cd pysisyphus/nix
    # Executing nix-build for the first time may take a while as your local Nix store
    # will be populated.
    # Please grab a coffee while the command runs or start your own coffee plantation
    # and return after you picked your first crop of coffee berries ☕.
    nix-build

You will likely depend on closed source software (ORCA, Turbomole, Gaussian, ...) , which is not freely redistributable. If the binary/source archives of these programs are missing from the Nix-store, the installation process will interrupt and tell you how to provide the required files. So it's a good idea to investigate the output of `nix-build` from time to time check, if manual intervention is required.

Running pysisyphus with Nix
===========================

You can now make pysisyphus available to your user environment (not recommended) by

.. code-block:: bash

    nix-env -f default.nix -i

or launch a interactive `nix-shell`_ with pysisyphus by

.. code-block:: bash

   nix-shell

from within the pysisyphus repository.

Do not be confused if the commands of the underlying quantum chemistry codes are not available. They are made available to directly to the pysisyphus entry point, but not necessarily to your shell.

.. _`Nix package manager`: https://nixos.org/download.html
.. _`NixOS-QChem`: https://github.com/markuskowa/NixOS-QChem
.. _`nix-shell`: https://nixos.org/nix/manual/#sec-nix-shell
.. _`nix manual`: https://nixos.org/manual/nix/stable/
.. _`Nix Pills`: https://nixos.org/guides/nix-pills/index.html
