Worked example
==============

This page provides a general overview of pysisyphus's capabilities for optimizing ground
state reaction paths for the tris-pericyclic reaction between a tropone derivative and
dimethylfulvene. You can find a short discussion of the reaction on Steven M. Bacharachs
`blog`_. Instead of using DFT as in the original `publication`_ we'll employ the semi-emperical
tight binding code `xtb`_ as it allows us to run all necessary calculations on common
desktop hardware in little time. Whereas the paper discusses the importance of
dynamical effects on the reaction outcome we'll neglect them here for now.

.. _`blog`: http://comporgchem.com/blog/?p=4329
.. _`publication`: https://pubs.acs.org/doi/10.1021/jacs.8b12674
.. _`xtb`: https://github.com/grimme-lab/xtb

.. _reaction_fig:
.. figure:: /images/pericyclic_reactants.png
    :width: 1000
    :alt: Reactants of tris-pericyclig reaction.

    Educts, TS and products of the tris-pericyclic reaction between a tropone derivative and
    dimethylfulvene.

To model this reaction we have to create educt and product geometries using the molecular
editor of our choice, e.g., `avogadro`_ or `TmoleX`_. Given a reasonable educt geometry a
suitable product geometry can be constructed by decreasing the distance between the
two ring systems and optimizing the result directly by `xtb`_. Sometimes it may be easier
to start with the product geometry, from which the educt can be obtained.

.. _`avogadro`: https://avogadro.cc/
.. _`TmoleX`: https://www.3ds.com/products-services/biovia/products/molecular-modeling-simulation/solvation-chemistry/turbomoler/

Given educts (`min_xtbopt.xyz`) and product (`prod_xtbopt.xyz`) we can create our YAML input.
Our goal is to obtain the barrier heights for the reaction shown in :numref:`reaction_fig` by
means of growing string (GS) calculation, transition state (TS) optimization, intrinsic reaction coordinate (IRC) calculation and subsequent optimization of the IRC endpoints. The full input is shown below. We'll go over it step by step.

.. code-block:: yaml

    geom:
     type: dlc
     fn: [min_xtbopt.xyz, prod_xtbopt.xyz]
    preopt:
    cos:
     type: gs
     max_nodes: 18
    opt:
     type: string
     align: False
     max_cycles: 20
    tsopt:
     type: rsirfo
     do_hess: True
     hessian_recalc: 5
    irc:
     type: eulerpc
     rms_grad_thresh: 0.0005
    endopt:
    calc:
     type: xtb
     pal: 6

The desired coordinate system and the file names of the input geometries are given in the
:code:`geom` block.

.. code-block:: yaml

    geom:
     type: dlc
     fn: [min_xtbopt.xyz, prod_xtbopt.xyz]

Here we want to run the GS in delocalized internal coordiantes (DLC), which is the
preferred way. Instead of :code:`type: dlc` Cartesian coordinates could be used by
:code:`type: cart`. If DLCs fail, Cartesian coordinates should be tried.


.. code-block:: yaml

    preopt:

The :code:`preopt` block is given without any additional keywords, so sane defaults will
be used for the preoptimization of educt and product geometries (Rational function
optimization (RFO) in redundant internal coordinates). In our case preoptimization is not
strictly necessary, as we already preoptimized the geometries using `xtb`_. But in general
if is advised to span chain of states (COS) like GS or nudged elastic band (NEB) between
stationary points on the potential energy surface.
If the educts and products are NOT covalently bound it may be a good idea to restrict the
number of preoptimization cycles to a small number (:code:`max_cycles: 10`), as these
optimizations are sometimes hard to converge. Please see :ref:`Optimization of Minima`
for a list of possible keywords in the :code:`preopt` block.

.. code-block:: yaml

    cos:
     type: gs
     max_nodes: 18

The :code:`cos` block configures COS calculations. Here we request a GS (:code:`gs`)
with 18 nodes (images) between educt and product for a total of 20 nodes. By default
the first and last  node (educt and product) are kept fixed throughout the optimization.
Please see :ref:`Chain Of States Methods` for further information on COS methods.

.. code-block:: yaml

    opt:
     type: string
     align: False
     max_cycles: 20

COS/GS optimization is controlled via the :code:`opt` block. For GS one should always
use :code:`type: string`. In internal coordinates we disable automated alignment of geometries,
as it is not needed. We also restrict the number of optimization cycles to 20. Default is 50.
The chosen optimizer will do steepest descent (SD) steps when the string grew in the previous
cycle, otherwise conjugate gradient (CG) steps are used. When the GS is fully grown/connected
the optimizer will use limited-memory Broyden-FletcherGoldfarb-Shanno (L-BFGS).

.. code-block:: yaml

    tsopt:
     type: rsirfo
     do_hess: True
     hessian_recalc: 5

When the GS is converged its the energy image (HEI) is determined by cubic splining
and used as guess for a classical TS optimization using restricted step image RFO (RSIRFO)
:code:`do_hess: True` requests a frequency calculation after the TS optimization.
The Hessian is recalculated every 5th step. When the Hessian for the chosen computational
method is reasonably cheap it is a good idea to recalculate it periodically.
Between recalculations it is updated using the Bofill-update.

.. code-block:: yaml

    irc:
     type: eulerpc
     rms_grad_thresh: 0.0005

The :code:`irc` blocks controls the IRC integration. By default the Euler-predictor-corrector
(EulerPC) integrator is used. Integration is terminated when the root-mean-square of the
gradient is equal to or less than 0.0005 au. Possible inputs are given in
:ref:`Intrinsic Reaction Coordinate (IRC)`.


.. code-block:: yaml

    endopt:

Similar to :code:`preopt` the :code:`endopt` will be executed with default arguments. It
is used to optimize the IRC endpoints to stationary points and enables printing of
additional information like RMS deviation of atomic positions (RMSD) between optimized
endpoints and initial geometries.

.. code-block:: yaml

    calc:
     type: xtb
     pal: 6

The :code:`calc` block configures the level of theory used in energy/gradient/Hessian
calculations. Here we chose `xtb`_ and requested 6 CPU cores. Additional inputs for
xtb can be found in the :ref:`xtb module documentation <pysisyphus.calculators.XTB module>`

With everything set up we are ready to actually execute pysisyphus! Assuming the above
YAML is saved to `01_pericyclic.yaml` just run

.. code-block:: bash

    pysis 01_pericyclic.yaml | tee pysis.log

By default pysisyphus prints to STDOUT so you have to explicitely capture STDOUT. We use
:code:`tee` so everything is logged to a file and printed simulatenously.
A lot of files and output will be produced so we will go over everything slowly.

.. code-block::

                               d8b                            888
                               Y8P                            888
                                                              888
    88888b.  888  888 .d8888b  888 .d8888b  888  888 88888b.  88888b.  888  888 .d8888b
    888 "88b 888  888 88K      888 88K      888  888 888 "88b 888 "88b 888  888 88K
    888  888 888  888 "Y8888b. 888 "Y8888b. 888  888 888  888 888  888 888  888 "Y8888b.
    888 d88P Y88b 888      X88 888      X88 Y88b 888 888 d88P 888  888 Y88b 888      X88
    88888P"   "Y88888  88888P' 888  88888P'  "Y88888 88888P"  888  888  "Y88888  88888P'
    888           888                            888 888
    888      Y8b d88P                       Y8b d88P 888
    888       "Y88P"                         "Y88P"  888                            

    Version 0.4.3.post1+443.ga05d1aa4 (Python 3.8.5, NumPy 1.19.2, SciPy 1.5.2)
    Git commit a05d1aa4f63f480f4d61e2a15cf3a5b455097730
    Executed at Mon Oct  5 15:22:34 2020 on 'your_hostname'

    If pysisyphus benefitted your research please cite:

        https://doi.org/10.1002/qua.26390

    Good luck!

You will be greeted by a banner and some information about your current installation,
which hopefully aids in reproducing your results later on, if needed.
