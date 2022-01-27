let
  pkgs = import ./pkgs.nix;
  pysis = { pysisyphus = import ./default.nix { inherit fullTest; }; };
  allPkgs = pkgs // pysis;
in with allPkgs; mkShell {
  buildInputs = [ python3.pkgs.qcengine qchem.gamess-us ];
}
