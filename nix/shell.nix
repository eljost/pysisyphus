let pkgs = import ./pkgs.nix;
    pysisyphus = pkgs.python3.pkgs.pysisyphus.overrideAttrs (_: {
      doCheck = false;
      doInstallCheck = false;
    });
in with pkgs; mkShell {
  buildInputs = [pysisyphus ]
    ++ pysisyphus.nativeBuildInputs
    ++ pysisyphus.buildInputs
    ++ pysisyphus.propagatedBuildInputs
  ;
}
