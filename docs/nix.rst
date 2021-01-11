Installation via Nix
********************

An easy and reliable way to get a complete pysisyphus installation, including additional quantum chemistry (QC) codes, is possible through the `Nix package manager`_. Required Nix files for the QC codes and their dependencies (linear algebra & MPI libraries etc.) are provided by the NixWithChemistry_ overlay.

Nix installations are fully reproducible and extremely simple to accomplish for any potential pysisyphus user, as most of the necessary configuration is already described in the respective Nix files.

The `Nix Pills`_ series provides an informative introduction to the general concepts behind Nix.

Prerequisites
=============

Please follow the installation instructions in the `nix manual`_ and additional hints from Chemix to set up a working Nix installation on your system. You'll probably have to allow unfree packages.

.. code-block:: bash
    nix-channel --add https://nixos.org/channels/nixos-20.09 nixos
    nix-channel --update
    export NIXPKGS_ALLOW_UNFREE=1

Binary Caching
--------------

Binary caching of Nix dependencies heavily speeds up build times by getting the heavy dependencies (like Psi4) prebuilt from the cache (usually hours). A binary cache for the open source components for a common workstation configuration is provided:

- SSE exentensions up to SSE4.2
- AVX extensions up to AVX2
- BLAS and LAPACK is provided by MKL
- MVAPICH2 is used as the common MPI provider
- The network for MPI parallelisation is ethernet and the loopback interface

The caching is provided by the Cachix_ service and needs to be enabled for your nix installation:

.. code-block:: bash
    nix-env -iA nixos.cachix
    cachix use chemix

Cachix will (necessarily) add a new :code:`substituters` and :code:`trusted-public-keys` to your :code:`$HOME/.config/nix/nix.conf` if executed as normal user or :code:`/etc/nix/nix.conf` if executed as root.

If you intend to run `nix-build` as normal user, this user *must be* added to `trusted-users` in :code:`/etc/nix/nix.conf` or the caching settings are *ignored*
and everything will be built on your computer.

Configuration
=============

Only minor configuration is required for a successful pysisyphus installation with Nix: Clone the pysisyphus repository from github, initialize the NixWithChemistry_ submodule and adapt `nix/nixwithchemistry/config.nix` as necessary.

.. code-block:: bash

    git clone https://github.com/eljost/pysisyphus.git $install_dir
    cd $install_dir
    git submodule init nix/nixwithchemistry
    git submodule update
    # Adapt the configuration file with the $EDITOR of your choice (nano/vim/kate/gedit...)
    edit nix/nixwithchemistry/config.nix

Please adapt the following lines in :code:`config.nix`:
    - The values for proprietary codes, e.g. :code:`gaussian`, :code:`mrcc`, :code:`orca`, :code:`turbomole` ..., can be set to :code:`true` if you have access to their installers and plan to use them.
    - Adapt the available CPU extensions on your machine (you can list them with :code:`grep "flags" /proc/cpuinfo`)
    - :code:`network` should be adapted if you intend to run pysisyphus on a HPC cluster. Possible values are: :code:`network = "ethernet"`, :code:`network = "omnipath"` or :code:`network = "infiniband"`. If you are running pysisyphus on local machine keep :code:`network = "ethernet"`.
    - :code:`mpi` is set to :code:`"mvapich2"` by default. There  is usually no reason to change this but you could also use :code:`"openmpi"`. This does **not** mean that the selected MPI version will be used in all programs. MRCC requires IntelMPI and GAMESS-US OpenMPI, no matter what you select here.
    - :code:`blas = "mkl"` can be updated to :code:`blas = "openblas"` if you have an AMD CPU, but it is probably better to set the :code:`MKL_DEBUG_CPU_TYPE=5` environment variable.

While enabling the proprietary components in :code:`config.nix` usually does not have an impact on binary caching, changing the CPU optimisation flags, the MPI implementation or network settings will usually trigger more rebuilds as caching is not setup for those configurations.

Installation
============

After you updated :code:`nix/nixwithchemistry/config.nix`, you can start building/fetching the actual QC codes and all dependencies.

The initial `nix-build` took about 3.5 h on my laptop, equipped with a i7-8750H (6 physical cores), 16 GB RAM and a NVMe SSD and requires about 15 GB disk space in the local Nix store (only ORCA was enabled, Gaussian and Turbomole were disabled).

.. code-block:: bash

    cd $install_dir/nix
    # Executing nix-build for the first time may take a while as your local Nix store
    # will be populated.
    # Please grab a coffee while the command runs or start your own coffee plantation
    # and return after you picked your first crop of coffee berries ☕.
    nix-build

You will likely depend on closed source software (ORCA, Turbomole, Gaussian, ...) , which is not freely redistributable. If the binary/source archives of these programs are missing from the Nix-store, the installation process will interrupt and tell you how to provide the required files. So it's a good idea to investigate the output of `nix-build` from time to time check, if manual intervention is required.

Running pysisyphus with Nix
===========================

You can now make pysisyphus available to your user environment by

.. code-block:: bash

    nix-env -f default.nix -i

or launch a `nix-shell`_ with pysisyphus by

.. code-block:: bash

   nix-shell --pure

or use :code:`nix run`

.. code-block:: bash

    nix run

**WARNING** In case of :code:`nix run` the resulting shell will not be pure. Depending on your system configuration conda/pip/... packages and configurations from the system might leak in. You are definitely safe with :code:`nix-shell --pure`.

Do not be confused if the commands of the underlying quantum chemistry codes are not available. They are made available to directly to the pysisyphus entry point, but not necessarily to your shell.

.. _`Nix package manager`: https://nixos.org/download.html
.. _NixWithChemistry: https://gitlab.com/theoretical-chemistry-jena/nixwithchemistry
.. _`nix-shell`: https://nixos.org/nix/manual/#sec-nix-shell
.. _`nix manual`: https://nixos.org/manual/nix/stable/
.. _`Nix Pills`: https://nixos.org/guides/nix-pills/index.html
.. _Cachix: https://cachix.org/
