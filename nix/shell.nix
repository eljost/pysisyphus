let
  sources = import ./sources.nix;
  qchem = import sources.NixOS-QChem;
  nixpkgs = import sources.nixpkgs {
    overlays = [ qchem ];
    allowUnfree = true;
  };
  pysis = { pysisyphus = import ./default.nix; };
  allPkgs = nixpkgs // pysis;
in
  with allPkgs;
  mkShell {
    buildInputs = [ pysisyphus ];
  }
