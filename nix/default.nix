{ fullTest ? false
, postOverlays ? []
} :

let pkgs = import ./pkgs.nix { inherit postOverlays; };
in with pkgs; qchem.python3.pkgs.callPackage ./pysisyphus.nix {
  orca = qchem.orca; # Uses the screenreader if not given

  # Perform full testing?
  inherit fullTest;
}
