let pkgs = import ./pkgs.nix;

in with pkgs; qchem.python3.pkgs.callPackage ./pysisyphus.nix {
  # Enable quantum chemistry codes by setting them to true.
  # The default settings in ./pysisyphus.nix enable only free codes
  #enableJmol = true;
  #enableMultiwfn = true;
  #enableXtb = true;
  #enableOpenmolcas = true;
  #enablePsi4 = true;
  #enableWfoverlap = true;
  #enableNwchem = true;
  #enableOrca = true;
  #enableTurbomole = true;
  #enableGaussian = true;
  #enableCfour = true;
  #enableMolpro = true;
  #enableGamess = true;

  # Usually no need to change.
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
  gamess-us = qchem.gamess-us.override { enableMpi = false; };
}
