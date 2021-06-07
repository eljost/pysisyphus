{ fullTest ? false
, postOverlays ? []
} :

let pkgs = import ./pkgs.nix { inherit postOverlays; };
in with pkgs; qchem.python3.pkgs.callPackage ./pysisyphus.nix {
  multiwfn = qchem.multiwfn;
  xtb = qchem.xtb;
  wfoverlap = qchem.wfoverlap;
  nwchem = qchem.nwchem;

  # Perform full testing?
  inherit fullTest;

  # Uncomment below to enable optional engines.
  orca = qchem.orca;
  turbomole = qchem.turbomole;
  cfour = qchem.cfour;
  molpro = qchem.molpro;
  gaussian = qchem.gaussian;
}
