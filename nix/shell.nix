{ fullTest ? false
, postOverlays ? []
} :

let
  pkgs = import ./pkgs.nix { inherit postOverlays; };
  pysis = { pysisyphus = import ./default.nix { inherit fullTest; }; };
  allPkgs = pkgs // pysis;
in
  with allPkgs;
  mkShell {
    buildInputs = [ pysisyphus ];
  }
