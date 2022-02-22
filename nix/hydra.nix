let pkgs = (import ./pkgs.nix) // { config.allowUnfree = true; };
in
{ pysisyphus = pkgs.python3.pkgs.pysisyphus.override {
    enableMultiwfn = true;
    enableOrca = true;
    enableTurbomole = true;
    enableCfour = true;
    enableGamess = true;
  };
}
