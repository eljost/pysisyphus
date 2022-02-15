self: super: {
  pysisyphus = super.python3.pkgs.toPythonApplication self.python3.pkgs.pysisyphus;

  python3 = super.python3.override { packageOverrides = pself: psuper: {
    pysisyphus = super.qchem.python3.pkgs.callPackage ./pysisyphus.nix {
      inherit (super.qchem)
        multiwfn
        xtb
        molcas
        psi4
        wfoverlap
        nwchem
        orca
        turbomole
        gaussian
        cfour
        molpro
        gamess-us
      ;
    };
  };};
}
