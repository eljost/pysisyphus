Minimization with Microiterations
*********************************

Pysisyphus allows layered optimizations using  microiterations via the
`LayerOpt` class. Required energies and its derivatives (gradients, Hessians)
can either be calculated with pysisyphus' own ONIOM implementation or supplied
through Unix sockets via the i-PI protocol.

In ONIOM calculations a small (model) system is treated using a high level of theory,
while its surroundings (real system) are described at a lower level of theory. There
may also be several (optional) intermediate layers. Please see the figure below.

.. figure:: _static/oniom.svg
  :scale: 50%
  :width: 100%
  :align: center
  :alt: General layer structure in a ONIOM calculation.

  General layer structure in a ONIOM calculation. The indices of the respective layers
  are given as red numbers in the top right corner of each layer.

Searching stationary points of such systems using classical optimization approaches with
simultanous relaxation of model and real system can be computationally quite costly.
In cases where the model system converges before the real system, many unnecessary
model system energy and gradients evaluations may be required to fully relax
the total system, as one step requires gradient evaluations in all layers.

Optimization with microiterations offers a way to potentially reduce the number
of costly calculations in the innermost layer. By fully relaxing the outer layers
before taking a step in the innermost layer, superfluous gradient evaluations can be avoided.
As the overall calculation time of an ONIOM energy/gradient is usually dominated by the
innermost layer, great compuational saving can by realized by decreasing the
required number of such calculations. 

In the present approach, outer layers are optimized in Cartesian coordinates with an
economic optimizer like (preconditioned) limited memory BFGS (L-BFGS) or conjugate
gradient (CG).
As L-BFGS avoids operations with steep scaling costs like matrix diagonalization,
relaxation of large system bcomprising thousands of atoms becomes possible with
negligible computational overhead.
In contrast to the outer layers, the inner layer is usually optimized using internal
coordinates and an optimizer that utilizes an explicit Hessian matrix as this often
allows for the most efficient optimizations.
As energy and gradient evaluations in the outer layers are usually (very) cheap
the inferior performance of optimizations in Cartesian coordinates should not lead
to an overall runtime increase.

Input keywords
--------------

A layered optimization is requested via ``type: layer`` in the ``opt:`` section. When
pysisyphus' native ONIOM calculator is to be used the appropriate input has to be
given in the ``calc:`` section, e.g., layer composition and the different calculators.
If energies & gradients are sent via sockets and the i-PI-protocol ``calc:`` can be
left empty.

.. code:: yaml

   opt:
    type:
     layer:

Example using pysisyphus' ONIOM
-------------------------------

.. literalinclude :: ../examples/opt/21_xtb_layeropt/21_xtb_layeropt.yaml
   :language: yaml

Example with sockets & i-PI-protocol
------------------------------------

Whereas pysisyphus can figure out the layer composition when its own ONIOM calculator
is used, the user has to specify the layer structure when using sockets & the i-PI-protocol.

.. automodule:: pysisyphus.optimizers.LayerOpt
    :members:
    :undoc-members:
    :show-inheritance:
