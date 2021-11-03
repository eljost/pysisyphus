Excited State Optimization
**************************

pysisyphus offers excited state (ES) tracking, to follow diadiabatic states
along an optimization. ES tracking is enabled by putting `track: True` in the
`calc` section of your YAML input. In ES optimizations pysisyphus may use
additional programs (wfoverlap, Multiwfn, jmol) if requested. Please see the
:ref:`installation instructions<pysisrc-label>` for information on how to set
them up.

**Please consider reading the relevant sections (2 and 3; 4 discusses examples) of
the** `pysisyphus paper <https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26390>`_.
Additionally the user should think about the relevance between equilibrium/non-equilibrium
solvation when calculating ES-gradients with implicit solvation. Serveral programs,
e.g., Gaussian use equilibrium solvation when doing ES-optimization. It is in the
responsibility of the user to add the relevant keywords when using pysisyphus, e.g.,
:code:`EqSolv` for Gaussian.

YAML Example
------------

A bare-bone input for the S\ :sub:`1` optimization of the 1H-amino-keto
tautomer of cytosin at the TD-DFT/PBE0/def2-SVP level of theory using ORCA is
shown below. The full example is found
`here <https://github.com/eljost/pysisyphus/tree/master/examples/opt/06_orca_cytosin_s1_opt>`_.

.. code:: yaml

    opt:
    calc:
     type: orca
     keywords: pbe0 def2-svp rijcosx def2/J def2-svp/C
     # Calculate 2 ES by TD-DFT, follow the first one
     blocks: "%tddft nroots 2 iroot 1 tda false end"
     charge: 0
     mult: 1
     pal: 4
     mem: 2000
     # ES-tracking related keywords follow from here
     # Enable ES-tracking, this is important.
     track: True
     # Track ES by transition density overlaps
     ovlp_type: tden
    geom:
     type: redund
     fn: cytosin.xyz

Additional keywords are possible in the `calc` section. The default values are shown
below.

.. code:: yaml

    calc:
     # Controls calculation of charge-density-differences cubes and rendering
     # Cubes are calcualted by Multiwfn, rendering is handled by jmol.
     #
     # Possible values are: (None, calc, render).
     #  None: No cube calculation/rendering
     #  calc: Cube calculation by Multiwfn.
     #  render: Same as 'calc' and cubes are then rendered by jmol.
     cdds: None
     # Overlap type. Using 'wf' requires the external wfoverlap binary. The remaining
     # options are implemented directly in pysisyphus.
     #
     # Possible values are (wf, tden, nto, nto_org). 
     #  wf: Wavefunction overlaps using the external wfoverlap program.
     #  tden: Transition density matrix overlaps.
     #  nto: Natural transition orbital overlaps.
     #  nto_org: Natural transition orbital overlaps as described by García.
     ovlp_type: wf
     # Controls the reference cycle that is used in the overlap calculation. The default
     # 'adapt' is recommended.
     #
     # Possible values are (first, previous, adapt)
     #  first: Keep first calculation as reference cycle. Reliable when only minor
     #         geometrical changes are expected.
     #  previous: Use previous cycle as reference.
     #  adapat: Use adaptive algorithm. Please see the pysisyphus paper for a discussion
     #          at the end of section 3.
     ovlp_with: adapt
     # Thresholds controlling the update of the reference cycle.
     # The first number specifies the minimum overlap that must be exceeded, for an update
     # of the reference cycle. Assuming a value of 0.5 (50 %), the reference cycle update
     # is skipped, if the overlaps between the current states and the reference state don't
     # exceed 50 %.
     # The last two numbers define an interval for the ratio between the second highest
     # overlap, and the highest overlap. If the ratio is small, e.g., below 0.3, then both
     # states are sufficiently different, and no reference cycle update is needed. If the ratio
     # is bigger (> 0.6), then the states are quite similar, and an update is currently not
     # advised.
     # Possible values [three positive floats between 0. and 1.]
     adapt_args: [0.5, 0.3, 0.6]
     # Explicitly calculate the AO-overlap matrix in a double molecule calculation. Only
     # supported by Turbomole and Gaussian calculators. If False, the approximate AO
     # overlap matrix is reconstructed from inverting the MO-coefficient matrix.
     #
     # Possible values: (True, False)
     double_mol: False
     # Absolute CI-coefficients below this threshold are ignored in the overlap calculation.
     #
     # Possible values: positive float
     conf_thresh: 0.0001
     #
     # nto/natural transition orbital specific
     #
     # Number of NTOs to consider in the overlap calculation. Only relevant for 'nto'
     # and 'nto_org' ovlp_types.
     #
     # Possible values: positive integer
     use_ntos: 4
     # Dynamically decide on number of NTOs according to their participation ratio. Only
     # relevant for 'nto_org'
     #
     # Possible values: boolean
     pr_nto: False
     # 
     # wfoverlap/wavefunction overlaps specific
     # 
     # Number of core orbitals to neglect in a wfoverlap calculation. Only relevant
     # for the 'wf' ovlp_type. Must be >= 0.
     #
     # Possible values: positive integer
     ncore: 0
     #
     # tden/transition density matrix specific
     #
     # Controls which set of MO coefficients (at current cycle, or the reference cycle)
     # is used to recover the AO overlap matrix.
     #
     # Possible values: (ref, cur)
     mos_ref: cur
     # Controls whether the set of MO coefficents that was NOT used for recovering the AO
     # overlap matrix is re-normalized, using the recovered AO overlap matrix. If set to
     # True and mos_ref = cur, then the MO coefficients at the reference cycle will be re-
     # normalized, and vice versa.
     #
     # Possible values: (True, False)
     mos_renorm: True

By brief reasoning it would seem that :code:`mos_ref: ref` and :code:`mos_renorm: True` are
more sensible choices, which is possibly true. Right now the present defaults are kept for
legacy reasons, and I'll update them after testing out the alternatives.

Please also see :ref:`Example - Excited State Tracking <Plotting ES optimizations>`
for possible visualizations when optimizing ES.

Optimization of Conical Intersections
-------------------------------------

pysisyphus implements the `projected gradient method using
an updated branching plane`_, as developed
by Maede, Ohno and Morokuma. Currently, CI-optimization is not enabled for YAML input.
An illustrative example is found in *tests/test_conic_intersect*.

.. _projected gradient method using an updated branching plane: https://pubs.acs.org/doi/pdf/10.1021/ct1000268
