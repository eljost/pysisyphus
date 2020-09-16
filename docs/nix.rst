Nix
***

The easiest way to get a complete pysisyphus installation with quantum chemistry dependencies is to use the `nix package manager`_.
A nix installation of pysisyphus comes with "batteries included" and contains the required quantum chemistry programs, provided by the NixWithChemistry_ overlay.
A nix installation avoids the problems of global installations and virtual environments, is portable and fully reproducible.

Prerequisites
=============

Nix needs to be installed on your system and working. Follow the installation instructions in the `nix manual`_ and possibly the additional hints from NixWithChemistry_.

Configuration
=============

Quantum chemistry software is a performance critical part of pysisyphus, might depend on your hardware configuration and its availability might depend on licenses.
Therefore a small amount of configuration work is necessary for the software from the NixWithChemistry_ repository.

Clone the code from github, initialise the NixWithChemistry_ submodule and make changes to the configuration.

.. code-block:: bash

    git clone https://github.com/eljost/pysisyphus.git $install_dir
    cd $install_dir
    git submodule init nix/nixwithchemistry
    git submodule update
    # Make changes to the configuration file.
    edit nix/nixwithchemistry/config.nix

In :code:`config.nix` make sure to have the following set correctly:
    - :code:`gaussian`, :code:`mrcc`, :code:`orca` and :code:`turbomole` should be set to :code:`false` if no license and source/binary archive is available.
    - The CPU extensions available on your machine should be set correctly. You can check your machine with :code:`grep "flags" /proc/cpuinfo`.
    - If you are not using CUDA enabled hardware, set :code:`cuda = false;`.
    - :code:`network` should be changed if pysisyphus is going to run on a HPC cluster, where you could change :code:`network = "ethernet"` to :code:`network = "omnipath"` or :code:`cuda = "infiniband"`. If you are running on local machine keep :code:`cuda = "ethernet"`.
    - :code:`mpi` is set to :code:`"mvapich2"` by default. There  is usually no reason to change this but you could also use :code:`"openmpi"`. This does **not** mean that the selected MPI version will be used in all programs. MRCC requires IntelMPI and GAMESS-US OpenMPI, no matter what you select there.
    - :code:`blas` can be set to :code:`"mkl"` or :code:`"openblas"`. On Intel CPUs MKL is more efficient, otherwise you should use OpenBLAS.


Installation
============

After you made the changes to :code:`nix/nixwithchemistry/config.nix`, you are ready to build.

.. code-block:: bash

    cd $install_dir/nix
    nix-build

Likely you will depend on closed source software, which is not redistributable.
If the binary/source archives of these programs are missing from the nix-store, the installation process will interrupt and will tell you how to add the necessary files to the nix-store.

You can now make pysisyphus available to your user environment by

.. code-block:: bash

    nix-env -f default.nix -i

or launch a `nix-shell`_ with pysisyphus available by

.. code-block:: bash

   nix-shell --pure

or use :code:`nix run`

.. code-block:: bash

    nix run

.. _`nix package manager`: https://nixos.org/download.html
.. _NixWithChemistry: https://gitlab.com/theoretical-chemistry-jena/nixwithchemistry
.. _`nix-shell`: https://nixos.org/nix/manual/#sec-nix-shell
.. _`nix manual`: https://nixos.org/manual/nix/stable/
