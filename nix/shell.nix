let pkgs = import ./pkgs.nix;
in with pkgs; mkShell {
  buildInputs = [ python3.pkgs.pysisyphus ]
    ++ python3.pkgs.pysisyphus.nativeBuildInputs
    ++ python3.pkgs.pysisyphus.buildInputs
    ++ python3.pkgs.pysisyphus.propagatedBuildInputs
  ;
}
