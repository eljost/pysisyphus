{ fullTest ? false
, postOverlays ? []
} :

let pkgs = import ./pkgs.nix { inherit postOverlays; };
in with pkgs; qchem.python3.pkgs.callPackage ./pysisyphus.nix {
  multiwfn = qchem.multiwfn;
  xtb = qchem.xtb;
  openmolcas = qchem.molcas;
  psi4 = qchem.python3.pkgs.psi4;
  wfoverlap = qchem.wfoverlap;
  nwchem = qchem.nwchem;
  orca = qchem.orca;
  turbomole = qchem.turbomole;
  gaussian = qchem.gaussian;
  cfour = qchem.cfour;
  molpro = qchem.molpro;

  # Perform full testing?
  inherit fullTest;
}
