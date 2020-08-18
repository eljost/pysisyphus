Excited state Optimization
**************************

pysisyphus offers excited state (ES) tracking, to follow diadiabatic states
along an optimization. ES tracking is enabled by putting `track: True` in the
`calc` section of your YAML input. In ES optimizations pysisyphus may use
additional programs (wfoverlap, Multiwfn, jmol) if requested. Please see the
:ref:`installation instructions<pysisrc-label>` for information on how to set
them up.

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
    xyz: cytosin.xyz
    coord_type: redund

Additional keywords are possible in the `calc` section. The default values are shown
below.

.. code:: yaml

    calc:
     # Controls calculation of charge-density-differences cubes and rendering
     # Cubes are calcualted by Multiwfn, rendering is handled by jmol.
     #
     # Possible values are: (None, calc, render).
     cdds: None
     # Overlap type. Using 'wf' requires the external wfoverlap binary. The remaining
     # options are implemented directly in pysisyphus.
     #
     # Possible values are (wf, tden, nto, nto_org). 
     ovlp_type: wf
     # Controls the reference cycle that is used in the overlap calculation. The default
     # 'adapt' is recommended.
     #
     # Possible values are (first, previous, adapt)
     ovlp_with: adapt
     # Explicitly calculate the AO-overlap matrix in a double molecule calculation. Only
     # supported by Turbomole and Gaussian calculators.
     double_mol: False
     # CI-coefficients below this threshold are ignored in the overlap calculation.
     conf_thresh: 0.0001
     # Number of NTOs to consider in the overlap calculation. Only relevant for 'nto'
     # and 'nto_org' ovlp_types.
     use_ntos: 4
     # Number of core orbitals to neglect in a wfoverlap calculation. Only relevant
     # for the 'wf' ovlp_type. Must be >= 0.
     ncore: 0

Please also see :ref:`Link <es-plotting-label>` for possible plotting options for ES tracking
and optimizations.
