let nixpkgs = import ./nixpkgs.nix;
    config = import ./nixwithchemistry/config.nix;
in  with nixpkgs; python3Packages.callPackage ./pysisyphus.nix {
      orca = if config.orca then orca else null;
      turbomole = if config.turbomole then turbomole else null;
      gaussian = if config.gaussian then gaussian else null;
      mopac = null;
    }
