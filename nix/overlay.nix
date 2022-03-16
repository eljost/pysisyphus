final: prev: {
  pysisyphus =
    prev.python3.pkgs.toPythonApplication final.python3.pkgs.pysisyphus;

  python3 = prev.python3.override {
    packageOverrides = pfinal: pprev: {
      pysisyphus = prev.qchem.python3.pkgs.callPackage ./pysisyphus.nix {
        inherit (prev.qchem)
          multiwfn xtb molcas psi4 wfoverlap nwchem orca turbomole gaussian
          cfour molpro gamess-us;
      };
    };
  };
}
