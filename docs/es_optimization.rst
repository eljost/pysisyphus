Excited state Optimization
**************************

pysisyphus offers excited state (ES) tracking, to follow the correct
diadiabatic state even when ES crossings occur along an optimization.
ES tracking is enabled by the `track: True` in the `calc` section of
the YAML input. An example input for the optimization of the
S\ :sub`1` of the 1H-amino-keto tautomer of cytosin using ORCA is shown
below:

.. code:: yaml

    opt:
     type: rfo
     thresh: gau
    calc:
     type: orca
     keywords: pbe0 def2-svp rijcosx def2/J def2-svp/C
     # Calculate 2 ES by TD-DFT, follow the first one
     blocks: "%tddft nroots 2 iroot 1 tda false end"
     # Enable ES-tracking
     track: True
     # Track ES by transition density overlaps
     ovlp_type: tden
     charge: 0
     mult: 1
     pal: 4
     mem: 2000
     cdds: render
    xyz: cytosin.xyz
    coord_type: redund

To be continued.
