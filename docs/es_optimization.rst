Excited state Optimization
**************************

pysisyphus offers excited state (ES) tracking, to follow diadiabatic states
along an optimization. ES tracking is enabled by putting `track: True` in the
`calc` section of your YAML input. In ES optimizations pysisyphus may use
additional programs (wfoverlap, Multiwfn, jmol) if requested. Please see the
:ref:`installation instructions<pysisrc-label>` for information on how to set
them up.

**Please consider reading the relevant sections (2 and 3; 4 discusses examples) of
the** `pysisyphus paper <https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26390>`_.

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
     # Number of NTOs to consider in the overlap calculation. Only relevant for 'nto'
     # and 'nto_org' ovlp_types.
     #
     # Possible values: positive integer
     use_ntos: 4
     # Number of core orbitals to neglect in a wfoverlap calculation. Only relevant
     # for the 'wf' ovlp_type. Must be >= 0.
     #
     # Possible values: positive integer
     ncore: 0

Please also see :ref:`Link <es-plotting-label>` for possible plotting options for ES tracking
and optimizations.
