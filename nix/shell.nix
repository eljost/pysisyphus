let
  nixpkgs = import ./nixpkgs.nix;
  pysis = { pysisyphus = import ./default.nix; };
  allPkgs = nixpkgs // pysis;
in
  with allPkgs;
  mkShell {
    buildInputs = [
      pysisyphus
    ];
  }
